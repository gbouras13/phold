"""
Sub-database (ACR / VFDB / CARD / NetFlax / DefenseFinder) output writers.

Each sub-DB filters merged_df for its phrog tag, joins with its metadata
table, drops the bookkeeping columns, and writes a TSV. When no hits
exist for a sub-DB, an empty placeholder file is touched.

Pure polars — output formatting (``write_csv(separator="\\t")``) is
byte-identical to the previous ``pd.DataFrame.to_csv(index=False,
sep="\\t")`` implementation; verified against the snapshot suite in
tests/unit/test_sub_db_outputs.py for all 5 sub-DBs.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import polars as pl
from loguru import logger

from phold.utils.util import touch_file


# Bookkeeping columns produced upstream that should never make it into a
# sub-DB TSV. They originate from the weighted-bitscore expansion and the
# topfunction merge in phold.results.topfunction.
_COMMON_DROP_COLS: Tuple[str, ...] = (
    "phrog",
    "function",
    "product",
    "annotation_method",
    "function_with_highest_bitscore_proportion",
    "top_bitscore_proportion_not_unknown",
    "head_and_packaging_bitscore_proportion",
    "integration_and_excision_bitscore_proportion",
    "tail_bitscore_proportion",
    "moron_auxiliary_metabolic_gene_and_host_takeover_bitscore_proportion",
    "DNA_RNA_and_nucleotide_metabolism_bitscore_proportion",
    "connector_bitscore_proportion",
    "transcription_regulation_bitscore_proportion",
    "lysis_bitscore_proportion",
    "other_bitscore_proportion",
    "unknown_function_bitscore_proportion",
)


@dataclass(frozen=True)
class _SubDB:
    """One sub-database's I/O spec."""
    phrog: str                  # value of the `phrog` column to filter on
    join_col: str               # column to rename `tophit_protein` to AND join on
    metadata_file: str          # filename under the phold database dir
    metadata_sep: str           # CSV/TSV separator for the metadata file
    drop_join_col: bool         # whether to drop `join_col` from the output too
    out_name: str               # output TSV name under sub_db_tophits/

    def output_path(self, output: Path) -> Path:
        return Path(output) / "sub_db_tophits" / self.out_name


_SUBDBS: Tuple[_SubDB, ...] = (
    _SubDB("acr",           "reference",         "acrs_plddt_over_70_metadata.tsv",          "\t", True,  "acr_cds_predictions.tsv"),
    _SubDB("vfdb",          "prot_id",           "vfdb_description_output.csv",              ",",  False, "vfdb_cds_predictions.tsv"),
    _SubDB("card",          "Protein Accession", "card_plddt_over_70_metadata.tsv",          "\t", False, "card_cds_predictions.tsv"),
    _SubDB("netflax",       "protein",           "netflax_annotation_table.tsv",             "\t", False, "netflax_cds_predictions.tsv"),
    _SubDB("defensefinder", "reference",         "defensefinder_plddt_over_70_metadata.tsv", "\t", True,  "defensefinder_cds_predictions.tsv"),
)


def _write_sub_db_output(merged: pl.DataFrame, spec: _SubDB, database: Path, output: Path) -> None:
    """Produce one sub-DB's TSV: filter -> join -> drop -> write, or touch_file if no hits."""
    out_path = spec.output_path(output)
    hits = merged.filter(pl.col("phrog") == spec.phrog).rename({"tophit_protein": spec.join_col})

    if hits.is_empty():
        touch_file(out_path)
        return

    metadata = pl.read_csv(
        Path(database) / spec.metadata_file,
        separator=spec.metadata_sep,
        infer_schema_length=10_000,
    )
    joined = hits.join(metadata, on=spec.join_col, how="left")

    drop = list(_COMMON_DROP_COLS) + ([spec.join_col] if spec.drop_join_col else [])
    # Only drop columns that exist (defensive — input frame may evolve)
    drop_existing = [c for c in drop if c in joined.columns]
    joined = joined.drop(drop_existing)

    joined.write_csv(out_path, separator="\t")


def create_sub_db_outputs(merged_df: pl.DataFrame, database: Path, output: Path) -> bool:
    """Write per-sub-database TSVs to ``output/sub_db_tophits/``.

    Args:
        merged_df: Per-CDS predictions table (the one populated by phold
                   compare/topfunction). Must contain at least a ``phrog``
                   column and a ``tophit_protein`` column.
        database:  Path to the phold database directory (contains the
                   per-sub-DB metadata tables).
        output:    Output root; ``output/sub_db_tophits/`` is created.

    Returns:
        True (always; failures raise rather than return False).
    """
    sub_db_tophits_dir = Path(output) / "sub_db_tophits"
    sub_db_tophits_dir.mkdir(parents=True, exist_ok=True)

    for spec in _SUBDBS:
        _write_sub_db_output(merged_df, spec, database, output)

    return True
