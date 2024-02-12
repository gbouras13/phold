"""
functions for creating vfdb, card, defensefinder and acr outputs
"""

from pathlib import Path

# imports
import pandas as pd
from loguru import logger

from phold.utils.util import touch_file


def create_sub_db_outputs(
    merged_df: pd.DataFrame, database: Path, output: Path
) -> bool:

    sub_db_tophits_dir = Path(output) / "sub_db_tophits"
    sub_db_tophits_dir.mkdir(parents=True, exist_ok=True)

    # acr df
    acr_df = merged_df[(merged_df["phrog"] == "acr")]
    acr_df = acr_df.rename(columns={"tophit_protein": "reference"})

    acr_metadata_path: Path = Path(database) / "acrs_plddt_over_70_metadata.tsv"
    acr_metadata_df = pd.read_csv(acr_metadata_path, sep="\t")

    columns_to_drop = [
        "phrog",
        "function",
        "product",
        "annotation_method",
        "function_with_highest_bitscore_percentage",
        "top_bitscore_percentage_not_unknown",
        "top_bitscore_percentage_not_unknown",
        "head_and_packaging_bitscore_percentage",
        "integration_and_excision_bitscore_percentage",
        "tail_bitscore_percentage",
        "moron_auxiliary_metabolic_gene_and_host_takeover_bitscore_percentage",
        "DNA_RNA_and_nucleotide_metabolism_bitscore_percentage",
        "connector_bitscore_percentage",
        "transcription_regulation_bitscore_percentage",
        "lysis_bitscore_percentage",
        "other_bitscore_percentage",
        "unknown_function_bitscore_percentage",
    ]

    acr_merged_output_path: Path = Path(sub_db_tophits_dir) / "acr_cds_predictions.tsv"

    #  if there is an acr hit
    if len(acr_df) > 0:
        acr_merged_df = pd.merge(acr_df, acr_metadata_df, on="reference", how="left")

        acr_merged_df = acr_merged_df.drop(
            columns=columns_to_drop + ["reference"],
        )
        acr_merged_df.to_csv(acr_merged_output_path, index=False, sep="\t")
    else:
        touch_file(acr_merged_output_path)

    # vfdb df
    vfdb_df = merged_df[(merged_df["phrog"] == "vfdb")]
    vfdb_df = vfdb_df.rename(columns={"tophit_protein": "prot_id"})

    vfdb_metadata_path: Path = Path(database) / "vfdb_description_output.csv"

    vfdb_metadata_df = pd.read_csv(vfdb_metadata_path, sep=",")
    vfdb_merged_output_path: Path = (
        Path(sub_db_tophits_dir) / "vfdb_cds_predictions.tsv"
    )

    # cleanup only if it has a hit
    if len(vfdb_df) > 0:
        vfdb_merged_df = pd.merge(vfdb_df, vfdb_metadata_df, on="prot_id", how="left")
        vfdb_merged_df = vfdb_merged_df.drop(columns=columns_to_drop)
        vfdb_merged_df.to_csv(vfdb_merged_output_path, index=False, sep="\t")

    else:
        touch_file(vfdb_merged_output_path)

    # card df
    card_df = merged_df[(merged_df["phrog"] == "card")]
    card_df = card_df.rename(columns={"tophit_protein": "Protein Accession"})

    card_metadata_path: Path = Path(database) / "card_plddt_over_70_metadata.tsv"
    card_metadata_df = pd.read_csv(card_metadata_path, sep="\t")
    card_merged_output_path: Path = (
        Path(sub_db_tophits_dir) / "card_cds_predictions.tsv"
    )

    # cleanup only if it has a hit
    if len(card_df) > 0:

        card_merged_df = pd.merge(
            card_df, card_metadata_df, on="Protein Accession", how="left"
        )
        card_merged_df = card_merged_df.drop(columns=columns_to_drop)
        card_merged_df.to_csv(card_merged_output_path, index=False, sep="\t")

    else:
        touch_file(card_merged_output_path)

    # defensefinder df
    defensefinder_df = merged_df[(merged_df["phrog"] == "defensefinder")]
    defensefinder_df = defensefinder_df.rename(columns={"tophit_protein": "reference"})

    defensefinder_metadata_path: Path = (
        Path(database) / "defensefinder_plddt_over_70_metadata.tsv"
    )
    defensefinder_metadata_df = pd.read_csv(defensefinder_metadata_path, sep="\t")
    defensefinder_merged_output_path: Path = (
        Path(sub_db_tophits_dir) / "defensefinder_cds_predictions.tsv"
    )

    # cleanup only if it has a hit
    if len(defensefinder_df) > 0:

        defensefinder_metadata_df = pd.merge(
            defensefinder_df, defensefinder_metadata_df, on="reference", how="left"
        )
        defensefinder_metadata_df = defensefinder_metadata_df.drop(
            columns=columns_to_drop + ["reference"]
        )
        defensefinder_metadata_df.to_csv(
            defensefinder_merged_output_path, index=False, sep="\t"
        )

    else:
        touch_file(defensefinder_merged_output_path)

    return True