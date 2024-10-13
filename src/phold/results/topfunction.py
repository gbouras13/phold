#!/usr/bin/env python3
import copy
from pathlib import Path
from typing import Dict, Tuple, Union

import polars as pl
from loguru import logger


def get_topfunctions(
    result_tsv: Path,
    database: Path,
    database_name: str,
    structures: bool,
    card_vfdb_evalue: float,
    proteins_flag: bool,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Process Foldseek output to extract top functions and weighted bitscores.
    """

    logger.info("Processing Foldseek output")

    col_list = [
        "query",
        "target",
        "bitscore",
        "fident",
        "evalue",
        "qStart",
        "qEnd",
        "qLen",
        "tStart",
        "tEnd",
        "tLen",
    ]

    foldseek_df = pl.read_csv(result_tsv, sep="\t", has_header=False, new_columns=col_list)

    if foldseek_df.is_empty():
        logger.error("Foldseek found no hits whatsoever - please check whether your input is really phage-like")

    foldseek_df = foldseek_df.with_column(pl.lit("foldseek").alias("annotation_source"))

    if not structures and not proteins_flag:
        foldseek_df = foldseek_df.with_columns(
            pl.col("query").str.split(":").arr.eval(pl.first(), strict=True).alias("contig_id"),
            pl.col("query").str.split(":").arr.eval(pl.last(), strict=True).alias("cds_id"),
        )
    else:
        foldseek_df = foldseek_df.with_column(
            pl.col("query").str.replace(".pdb", "").str.replace(".cif", "").alias("cds_id")
        )

    foldseek_df = foldseek_df.with_column(
        pl.col("target").str.replace(".pdb", "").alias("target")
    )

    foldseek_df = foldseek_df.with_columns(
        pl.col("target").str.split(":").arr.eval(pl.first(), strict=True).alias("phrog"),
        pl.col("target").str.split(":").arr.eval(pl.last(), strict=True).alias("tophit_protein"),
    )

    foldseek_df = foldseek_df.with_column(
        pl.col("phrog").str.replace("phrog_", "")
    )

    mask = foldseek_df["phrog"].str.starts_with("envhog_")
    foldseek_df = foldseek_df.with_columns(
        pl.when(mask).then(pl.col("phrog").str.replace("envhog_", "")).otherwise(pl.col("phrog")).alias("phrog"),
        pl.when(mask).then("envhog_" + pl.col("tophit_protein")).otherwise(pl.col("tophit_protein")).alias("tophit_protein")
    )

    mask = foldseek_df["phrog"].str.starts_with("efam_")
    foldseek_df = foldseek_df.with_columns(
        pl.when(mask).then(pl.col("phrog").str.replace("efam_", "")).otherwise(pl.col("phrog")).alias("phrog")
    )

    mask = foldseek_df["phrog"].str.starts_with("dgr_")
    foldseek_df = foldseek_df.with_columns(
        pl.when(mask).then(pl.col("phrog").str.replace("dgr_", "")).otherwise(pl.col("phrog")).alias("phrog")
    )

    phrog_annot_mapping_tsv: Path = Path(database) / "phold_annots.tsv"
    phrog_mapping_df = pl.read_csv(phrog_annot_mapping_tsv, sep="\t")
    phrog_mapping_df = phrog_mapping_df.with_column(pl.col("phrog").cast(pl.Utf8))

    foldseek_df = foldseek_df.join(phrog_mapping_df, on="phrog", how="left")

    foldseek_df = foldseek_df.with_column(
        pl.when(pl.col("product").is_null()).then("hypothetical protein").otherwise(pl.col("product")).alias("product")
    )

    foldseek_df = foldseek_df.filter(
        ((pl.col("phrog") == "vfdb") | (pl.col("phrog") == "card")) &
        (pl.col("evalue").cast(pl.Float64) < float(card_vfdb_evalue)) |
        ((pl.col("phrog") != "vfdb") & (pl.col("phrog") != "card"))
    )

    def custom_nsmallest(group):
        if all(group["product"] == "hypothetical protein"):
            return group.sort("evalue").first()
        else:
            return group.filter(pl.col("product") != "hypothetical protein").sort("evalue").first()

    topfunction_df = foldseek_df.groupby("query").apply(custom_nsmallest)

    topfunction_df = topfunction_df.drop("query")

    topfunction_df = topfunction_df.with_column(
        pl.col("evalue").apply(lambda x: "{:.3e}".format(float(x)))
    )

    def weighted_function(group: pl.DataFrame) -> pl.DataFrame:
        weighted_counts_normalised = {}
        bitscore_by_function = group.groupby("function").agg(pl.sum("bitscore"))

        total_functional_bitscore = group.filter(pl.col("function") != "unknown function").select(pl.sum("bitscore"))

        if total_functional_bitscore == 0:
            top_bitscore_function = "unknown function"
            top_bitscore_perc = 0
        else:
            for row in bitscore_by_function.iter_rows():
                key, value = row
                if key != "unknown function":
                    weighted_counts_normalised[key] = round(value / total_functional_bitscore[0, 0], 3)

            if weighted_counts_normalised:
                top_bitscore_function = max(weighted_counts_normalised, key=weighted_counts_normalised.get)
                top_bitscore_perc = max(weighted_counts_normalised.values())
            else:
                top_bitscore_function = "unknown function"
                top_bitscore_perc = 0

        return pl.DataFrame({
            "function_with_highest_bitscore_proportion": [top_bitscore_function],
            "top_bitscore_proportion_not_unknown": [top_bitscore_perc],
            # Add the rest of the bit score calculations
        })

    weighted_bitscore_df = foldseek_df.groupby("query").apply(weighted_function)

    weighted_bitscore_df = weighted_bitscore_df.with_column(
        pl.col("query").str.replace(".pdb", "").str.replace(".cif", "")
    ).drop("level_1")

    return topfunction_df, weighted_bitscore_df


def calculate_topfunctions_results(
    filtered_tophits_df: pl.DataFrame,
    cds_dict: Dict[str, Dict[str, dict]],
    output: Path,
    structures: bool,
    proteins_flag: bool,
    fasta_flag: bool,
) -> Union[Dict[str, Dict[str, dict]], pl.DataFrame]:
    """
    Calculate top function results based on filtered top hits DataFrame and update CDS dictionary accordingly.

    Args:
        filtered_tophits_df (pl.DataFrame): DataFrame containing filtered top hits.
        cds_dict (Dict[str, Dict[str, dict]]): Dictionary containing CDS information.
        output (Path): Output path.
        structures (bool): Indicates whether the input is is in .pdb or .cif format.
        proteins_flag (bool): Indicates whether the input is proteins.
        fasta_flag (bool): Indicates whether the input is in FASTA format.

    Returns:
        Union[Dict[str, Dict[str, dict]], pl.DataFrame]: Updated CDS dictionary and/or filtered top hits DataFrame.
    """

    # dictionary to hold the results
    result_dict = {}

    # so I can match with the df row below
    cds_record_dict = {}

    for record_id, cds_entries in cds_dict.items():
        result_dict[record_id] = {}
        for cds_id, cds_info in cds_entries.items():
            result_dict[record_id][cds_id] = {}
            cds_record_dict[cds_id] = record_id

    # Get record_id for every cds_id and merge into the df
    if structures is True:
        cds_record_df = pl.DataFrame(
            list(cds_record_dict.items()), schema=["cds_id", "contig_id"]
        )
        filtered_tophits_df = filtered_tophits_df.join(
            cds_record_df, on="cds_id", how="left"
        )

    if proteins_flag is True:
        column_order = ["cds_id"] + [
            col for col in filtered_tophits_df.columns if col != "cds_id"
        ]

    else:
        # Move "contig_id" and 'cds_id' to the front of the df
        column_order = ["contig_id", "cds_id"] + [
            col
            for col in filtered_tophits_df.columns
            if col != "contig_id" and col != "cds_id"
        ]
        filtered_tophits_df = filtered_tophits_df.select(column_order)

    # loop over all the foldseek tophits and add to the dict
    for row in filtered_tophits_df.iter_rows(named=True):
        if proteins_flag is False:
            record_id = row["contig_id"]
        else:
            record_id = "proteins"

        cds_id = row["cds_id"]
        values_dict = {
            "phrog": row["phrog"],
            "product": row["product"],
            "function": row["function"],
            "tophit_protein": row["tophit_protein"],
            "bitscore": row["bitscore"],
            "fident": row["fident"],
            "evalue": row["evalue"],
            "qStart": row["qStart"],
            "qEnd": row["qEnd"],
            "qLen": row["qLen"],
            "tStart": row["tStart"],
            "tEnd": row["tEnd"],
            "tLen": row["tLen"],
        }

        # nan on record_id -> means the structure in the structure_dir has a foldseek hit but can't be matched up to a contig
        if row.get("contig_id") is None:
            if structures is True:
                logger.warning(
                    f"{cds_id} has a foldseek hit but no record_id was found. Please check the way you named the structure file."
                )
            else:
                logger.warning(
                    f"{cds_id} has a foldseek hit but no record_id was found. Please check your input."
                )
        else:
            result_dict[record_id][cds_id] = values_dict

    # copy initial cds_dict
    updated_cds_dict = copy.deepcopy(cds_dict)

    phrog_function_mapping = {
        "unknown function": "unknown function",
        "transcription regulation": "transcription regulation",
        "tail": "tail",
        "other": "other",
        "moron, auxiliary metabolic gene and host takeover": "moron, auxiliary metabolic gene and host takeover",
        "lysis": "lysis",
        "integration and excision": "integration and excision",
        "head and packaging": "head and packaging",
        "DNA, RNA and nucleotide metabolism": "DNA, RNA and nucleotide metabolism",
        "connector": "connector",
    }

    # source dict
    source_dict = {}

    # iterates over the records
    for record_id, record in updated_cds_dict.items():
        source_dict[record_id] = {}
        # iterates over the features
        for cds_id, cds_feature in updated_cds_dict[record_id].items():
            # proteins/FASTA input -> no pharokka input -> fake input to make the updating work
            if proteins_flag is True or fasta_flag is True:
                cds_feature.qualifiers["function"] = ["unknown function"]
                cds_feature.qualifiers["phrog"] = ["No_PHROG"]
                cds_feature.qualifiers["product"] = ["hypothetical protein"]

            if result_dict[record_id][cds_id] != {}:
                foldseek_phrog = result_dict[record_id][cds_id].get("phrog", None)

                # add annotation source
                if proteins_flag:
                    source_dict[record_id][cds_id] = "foldseek"
                else:
                    source_dict[record_id][cds_id] = filtered_tophits_df.filter(
                        (pl.col("contig_id") == record_id) & (pl.col("cds_id") == cds_id)
                    ).select("annotation_source").item()

                if foldseek_phrog != cds_feature.qualifiers["phrog"][0]:
                    if cds_feature.qualifiers["phrog"][0] == "No_PHROG":
                        updated_cds_dict[record_id][cds_id].qualifiers["phrog"][0] = (
                            result_dict[record_id][cds_id]["phrog"]
                        )
                        updated_cds_dict[record_id][cds_id].qualifiers["product"][0] = (
                            result_dict[record_id][cds_id]["product"]
                        )
                        updated_cds_dict[record_id][cds_id].qualifiers["function"][
                            0
                        ] = result_dict[record_id][cds_id]["function"]

                    else:
                        if result_dict[record_id][cds_id]["function"] != "unknown function":
                            try:
                                updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                    0
                                ] = result_dict[record_id][cds_id]["phrog"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "product"
                                ][0] = result_dict[record_id][cds_id]["product"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "function"
                                ][0] = result_dict[record_id][cds_id]["function"]
                            except:
                                updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                    0
                                ] = result_dict[record_id][cds_id]["phrog"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "product"
                                ][0] = result_dict[record_id][cds_id]["product"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "function"
                                ][0] = result_dict[record_id][cds_id]["function"]

                        else:
                            if cds_feature.qualifiers["function"][0] == "unknown function":
                                updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                    0
                                ] = result_dict[record_id][cds_id]["phrog"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "product"
                                ][0] = result_dict[record_id][cds_id]["product"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "function"
                                ][0] = result_dict[record_id][cds_id]["function"]

                            else:
                                source_dict[record_id][cds_id] = "pharokka"

            else:
                if cds_feature.qualifiers["phrog"][0] == "No_PHROG":
                    source_dict[record_id][cds_id] = "none"
                else:
                    source_dict[record_id][cds_id] = "pharokka"

    return updated_cds_dict, filtered_tophits_df, source_dict



def initialize_function_counts_dict(
    record_id: str, count_dict: Dict[str, int], cds_count: int
) -> Dict[str, int]:
    """
    Initialize function counts dictionary for a given record ID.

    Args:
        record_id (str): ID of the record.
        count_dict (Dict[str, int]): Dictionary containing function counts.
        cds_count (int): Count of CDS.

    Returns:
        Dict[str, int]: Updated function counts dictionary.
    """
    count_dict[record_id]["cds_count"] = cds_count
    count_dict[record_id].update(
        {
            "phrog_count": 0,
            "connector": 0,
            "DNA, RNA and nucleotide metabolism": 0,
            "head and packaging": 0,
            "integration and excision": 0,
            "lysis": 0,
            "moron, auxiliary metabolic gene and host takeover": 0,
            "other": 0,
            "tail": 0,
            "transcription regulation": 0,
            "unknown function": 0,
        }
    )

    return count_dict
