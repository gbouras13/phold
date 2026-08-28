"""
NCBI feature-table (``.tbl``) writer — GenBank submission output.

Mirrors the ``.tbl`` that Pharokka emits (``pharokka.post_processing.create_tbl``)
so a Phold-reannotated genome can be handed to ``table2asn`` for GenBank
submission without a format conversion step (issue #137).

Format (tab-separated, per the NCBI feature-table spec):

    >Feature <contig_id>
    <start>\t<stop>\t<feature_key>
    \t\t\t<qualifier>\t<value>

``start``/``stop`` are 1-based and inclusive, and are written in *transcription*
order: for a minus-strand feature the larger coordinate comes first. That
ordering is how the feature table encodes strand — there is no separate strand
column — so it is the one piece of the format that must not be "tidied".
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl
from loguru import logger

from Bio.SeqFeature import SeqFeature

from phold.utils.util import atomic_write_path


# Qualifiers emitted for a CDS, in the order NCBI's examples use. Values come
# from per_cds_df, which write_genbank() has already normalised to 1-based
# inclusive coordinates (see the issue #75/#77 comment there).
_CDS_QUALIFIER_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("product", "product"),
    ("function", "function"),
    ("annotation_method", "inference"),
    ("transl_table", "transl_table"),
)

# Phold/Pharokka bookkeeping qualifiers that must never reach a submission
# file: either they are internal identifiers, or table2asn rejects them.
# ``ID``/``locus_tag`` are dropped because NCBI assigns locus tags at
# submission time; ``score``/``source``/``phase`` are pipeline provenance.
_NON_CDS_QUALIFIER_BLOCKLIST = frozenset(
    {
        "ID",
        "locus_tag",
        "score",
        "source",
        "phase",
        "translation",
        "transl_table",
        "phrog",
        "top_hit",
        "isotype",
    }
)

# Pharokka writes the tRNA name to a ``/trna`` qualifier, which is not a valid
# NCBI qualifier — ``product`` is, and it is mandatory for the tRNA feature key.
# ``write_genbank`` renames it in place, but only for records it walks, so do
# the same mapping here rather than depending on that call having happened
# first: a tRNA reaching table2asn without a product is rejected.
_QUALIFIER_ALIASES: Tuple[Tuple[str, str], ...] = (("trna", "product"),)

# Preferred emission order for the non-CDS qualifiers Pharokka carries. Anything
# not listed keeps dict order after these.
_NON_CDS_QUALIFIER_ORDER: Tuple[str, ...] = (
    "gene",
    "product",
    "inference",
    "ncRNA_class",
    "anticodon",
    "tag_peptide",
    "rpt_family",
    "rpt_type",
    "rpt_unit_range",
    "rpt_unit_seq",
    "db_xref",
    "note",
)


def _orient(start: int, end: int, strand: Optional[int]) -> Tuple[int, int]:
    """Return (start, stop) in transcription order.

    A minus-strand feature is written high-coordinate first; that swap *is*
    the strand annotation in a feature table.
    """
    if strand is not None and int(strand) == -1:
        return end, start
    return start, end


def _write_cds_features(handle, contig_df: pl.DataFrame) -> int:
    """Write every CDS row of one contig. Returns the number written."""
    written = 0
    for row in contig_df.iter_rows(named=True):
        start, stop = _orient(int(row["start"]), int(row["end"]), row.get("strand"))
        handle.write(f"{start}\t{stop}\tCDS\n")

        for column, qualifier in _CDS_QUALIFIER_COLUMNS:
            value = row.get(column)
            # A null product/function would emit a bare qualifier line that
            # table2asn treats as a parse error, so skip rather than write "".
            if value is None or str(value).strip() == "":
                continue
            handle.write(f"\t\t\t{qualifier}\t{value}\n")
        written += 1
    return written


def _non_cds_sort_key(item: Tuple[str, SeqFeature]) -> int:
    """Sort non-CDS features by start coordinate, tolerating odd locations."""
    _, feature = item
    try:
        return int(feature.location.start)
    except (AttributeError, TypeError):
        return 0


def _write_non_cds_features(
    handle, features: Dict[str, SeqFeature]
) -> int:
    """Write tRNA / tmRNA / CRISPR / ncRNA features carried over from the input.

    Phold does not predict these — they are passed through from the Pharokka
    (or NCBI/Bakta) input GenBank — but a submission-ready feature table has to
    carry them, otherwise the .tbl silently describes fewer features than the
    .gbk Phold wrote alongside it.
    """
    written = 0
    for _, feature in sorted(features.items(), key=_non_cds_sort_key):
        location = getattr(feature, "location", None)
        if location is None:
            continue

        # BioPython locations are 0-based half-open; the feature table is
        # 1-based inclusive, so the start needs +1 and the end is already
        # the inclusive last base.
        try:
            start = int(location.start) + 1
            end = int(location.end)
        except (TypeError, ValueError):
            logger.warning(
                f"Skipping feature with an unparseable location in the .tbl: {feature.type}"
            )
            continue

        start, stop = _orient(start, end, location.strand)
        handle.write(f"{start}\t{stop}\t{feature.type}\n")

        # Copy: the same SeqFeature objects are shared with the GenBank writer,
        # so the alias rewrite below must not mutate the caller's features.
        qualifiers = dict(feature.qualifiers or {})
        for source_key, target_key in _QUALIFIER_ALIASES:
            aliased = qualifiers.pop(source_key, None)
            if aliased is not None and target_key not in qualifiers:
                qualifiers[target_key] = aliased

        ordered = [k for k in _NON_CDS_QUALIFIER_ORDER if k in qualifiers]
        ordered += [k for k in qualifiers if k not in _NON_CDS_QUALIFIER_ORDER]

        for key in ordered:
            if key in _NON_CDS_QUALIFIER_BLOCKLIST:
                continue
            value = qualifiers[key]
            # BioPython stores qualifiers as lists when parsed from a GenBank
            # file but as bare strings when Phold built the feature itself.
            if isinstance(value, (list, tuple)):
                value = value[0] if value else None
            if value is None or str(value).strip() == "":
                continue
            handle.write(f"\t\t\t{key}\t{value}\n")
        written += 1
    return written


def write_tbl(
    per_cds_df: pl.DataFrame,
    non_cds_dict: Dict[str, Dict[str, SeqFeature]],
    contig_ids: List[str],
    prefix: str,
    output: Path,
) -> Path:
    """Write the NCBI feature table for every contig.

    Args:
        per_cds_df: The per-CDS dataframe returned by ``write_genbank`` —
            already 1-based inclusive, with contig_id / start / end / strand /
            product / function / annotation_method / transl_table.
        non_cds_dict: ``{contig_id: {feature_id: SeqFeature}}`` for the
            tRNA / tmRNA / CRISPR features passed through from the input.
        contig_ids: Contig order to emit, so the .tbl matches the .gbk.
        prefix: Output filename prefix.
        output: Output directory.

    Returns:
        Path to the written .tbl.

    Note:
        Every contig gets a ``>Feature`` header even when it has no features,
        because table2asn matches records by that header — omitting an empty
        contig makes it look absent from the submission rather than empty.
    """
    out_path: Path = Path(output) / f"{prefix}.tbl"

    has_cds = per_cds_df.height > 0 and "contig_id" in per_cds_df.columns
    total = 0

    with atomic_write_path(out_path) as tmp, open(tmp, "w") as f:
        for contig_id in contig_ids:
            f.write(f">Feature {contig_id}\n")

            if has_cds:
                contig_df = per_cds_df.filter(pl.col("contig_id") == contig_id)
                total += _write_cds_features(f, contig_df)

            total += _write_non_cds_features(f, non_cds_dict.get(contig_id, {}))

    logger.info(f"Wrote {total} features across {len(contig_ids)} contig(s) to {out_path}")
    return out_path
