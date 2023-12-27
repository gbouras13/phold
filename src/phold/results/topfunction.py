#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import copy
from loguru import logger


def get_topfunctions(
    result_tsv: Path, database: Path, database_name: str
) -> pd.DataFrame:
    logger.info("Processing Foldseek output.")
    col_list = [
        "query",
        "target",
        "foldseek_alnScore",
        "foldseek_seqIdentity",
        "foldseek_eVal",
        "qStart",
        "qEnd",
        "qLen",
        "tStart",
        "tEnd",
        "tLen",
    ]

    foldseek_df = pd.read_csv(
        result_tsv, delimiter="\t", index_col=False, names=col_list
    )

    foldseek_df[["record_id", "cds_id"]] = foldseek_df["query"].str.split(
        ":", expand=True, n=1
    )

    if database_name == "all_phrogs" or database_name == "all_envhogs":
        # split the first column
        foldseek_df[["phrog", "tophit_protein"]] = foldseek_df["target"].str.split(
            ":", expand=True, n=1
        )
    elif database_name == "all_phrogs_pdb" or database_name == "all_phrogs_reps":
        foldseek_df["phrog"] = foldseek_df["target"].str.replace(".pdb", "")
        foldseek_df["tophit_protein"] = None

    foldseek_df = foldseek_df.drop(columns=["target"])
    foldseek_df["phrog"] = foldseek_df["phrog"].str.replace("phrog_", "")
    foldseek_df["phrog"] = foldseek_df["phrog"].astype("str")

    # read in the mapping tsv
    phrog_annot_mapping_tsv: Path = Path(database) / "phrog_annot_v4.tsv"
    phrog_mapping_df = pd.read_csv(phrog_annot_mapping_tsv, sep="\t")
    phrog_mapping_df["phrog"] = phrog_mapping_df["phrog"].astype("str")

    # join the dfs
    foldseek_df = foldseek_df.merge(phrog_mapping_df, on="phrog", how="left")
    # Replace NaN values in the 'product' column with 'hypothetical protein'
    foldseek_df["product"] = foldseek_df["product"].fillna("hypothetical protein")
    foldseek_df["foldseek_eVal"] = foldseek_df["foldseek_eVal"].apply(
        lambda x: "{:.3e}".format(float(x))
    )

    # maybe dont need a sort
    def custom_nsmallest(group):
        # where all the
        if all(group["product"] == "hypothetical protein"):
            min_row_index = group["foldseek_eVal"].idxmin()
            # Get the entire row
            return group.loc[min_row_index]
        else:
            group = group[group["product"] != "hypothetical protein"]
            min_row_index = group["foldseek_eVal"].idxmin()
            # Get the entire row
            return group.loc[min_row_index]

    topfunction_df = (
        foldseek_df.groupby("query", group_keys=True)
        .apply(custom_nsmallest)
        .reset_index(drop=True)
    )

    # Remove the original 'query' column
    topfunction_df = topfunction_df.drop(columns=["query"])

    # scientific notation to 3dp
    topfunction_df["foldseek_eVal"] = topfunction_df["foldseek_eVal"].apply(
        lambda x: "{:.3e}".format(float(x))
    )

    return topfunction_df


def calculate_topfunctions_results(
    filtered_tophits_df: pd.DataFrame, cds_dict: dict, output: Path
) -> None:
    # Convert the DataFrame to a nested dictionary
    result_dict = {}

    # instantiate the unique contig ids
    unique_contig_ids = filtered_tophits_df["record_id"].unique()

    # loop over the contigs
    for record_id in unique_contig_ids:
        result_dict[record_id] = {}

    # loop over all the foldseek tophits
    for _, row in filtered_tophits_df.iterrows():
        record_id = row["record_id"]
        cds_id = row["cds_id"]
        values_dict = {
            "phrog": row["phrog"],
            "product": row["product"],
            "function": row["function"],
            "tophit_protein": row["tophit_protein"],
            "foldseek_alnScore": row["foldseek_alnScore"],
            "foldseek_seqIdentity": row["foldseek_seqIdentity"],
            "foldseek_eVal": row["foldseek_eVal"],
            "qStart": row["qStart"],
            "qEnd": row["qEnd"],
            "qLen": row["qLen"],
            "tStart": row["tStart"],
            "tEnd": row["tEnd"],
            "tLen": row["tLen"],
        }
        result_dict[record_id][cds_id] = values_dict

    # get counds
    # copy initial cds_dict
    updated_cds_dict = copy.deepcopy(cds_dict)

    original_functions_count_dict = {}
    new_functions_count_dict = {}
    combined_functions_count_dict = {}

    # iterates over the records
    for record_id, record in updated_cds_dict.items():
        # instantiate the functions dicts
        original_functions_count_dict[record_id] = {}
        new_functions_count_dict[record_id] = {}
        combined_functions_count_dict[record_id] = {}

        original_functions_count_dict[record_id]["cds_count"] = len(
            updated_cds_dict[record_id]
        )
        original_functions_count_dict[record_id]["phrog_count"] = 0
        original_functions_count_dict[record_id]["connector"] = 0
        original_functions_count_dict[record_id][
            "DNA, RNA and nucleotide metabolism"
        ] = 0
        original_functions_count_dict[record_id]["head and packaging"] = 0
        original_functions_count_dict[record_id]["integration and excision"] = 0
        original_functions_count_dict[record_id]["lysis"] = 0
        original_functions_count_dict[record_id][
            "moron, auxiliary metabolic gene and host takeover"
        ] = 0
        original_functions_count_dict[record_id]["other"] = 0
        original_functions_count_dict[record_id]["tail"] = 0
        original_functions_count_dict[record_id]["transcription regulation"] = 0
        original_functions_count_dict[record_id]["unknown function"] = 0

        new_functions_count_dict[record_id]["cds_count"] = len(
            updated_cds_dict[record_id]
        )
        new_functions_count_dict[record_id]["phrog_count"] = 0
        new_functions_count_dict[record_id]["connector"] = 0
        new_functions_count_dict[record_id]["DNA, RNA and nucleotide metabolism"] = 0
        new_functions_count_dict[record_id]["head and packaging"] = 0
        new_functions_count_dict[record_id]["integration and excision"] = 0
        new_functions_count_dict[record_id]["lysis"] = 0
        new_functions_count_dict[record_id][
            "moron, auxiliary metabolic gene and host takeover"
        ] = 0
        new_functions_count_dict[record_id]["other"] = 0
        new_functions_count_dict[record_id]["tail"] = 0
        new_functions_count_dict[record_id]["transcription regulation"] = 0
        new_functions_count_dict[record_id]["unknown function"] = 0
        new_functions_count_dict[record_id]["changed_phrogs"] = 0
        new_functions_count_dict[record_id]["same_phrogs"] = 0
        new_functions_count_dict[record_id]["different_phrogs"] = 0
        new_functions_count_dict[record_id]["foldseek_only_phrogs"] = 0
        new_functions_count_dict[record_id]["pharokka_only_phrogs"] = 0

        combined_functions_count_dict[record_id]["cds_count"] = len(
            updated_cds_dict[record_id]
        )
        combined_functions_count_dict[record_id]["phrog_count"] = 0
        combined_functions_count_dict[record_id]["connector"] = 0
        combined_functions_count_dict[record_id][
            "DNA, RNA and nucleotide metabolism"
        ] = 0
        combined_functions_count_dict[record_id]["head and packaging"] = 0
        combined_functions_count_dict[record_id]["integration and excision"] = 0
        combined_functions_count_dict[record_id]["lysis"] = 0
        combined_functions_count_dict[record_id][
            "moron, auxiliary metabolic gene and host takeover"
        ] = 0
        combined_functions_count_dict[record_id]["other"] = 0
        combined_functions_count_dict[record_id]["tail"] = 0
        combined_functions_count_dict[record_id]["transcription regulation"] = 0
        combined_functions_count_dict[record_id]["unknown function"] = 0

        # iterates over the features
        # maybe can add 3DI as a genbank feature eventually?
        for cds_id, cds_feature in updated_cds_dict[record_id].items():
            # if pharokka got a phrog
            if cds_feature.qualifiers["phrog"][0] != "No_PHROG":
                original_functions_count_dict[record_id]["phrog_count"] += 1
            # get original function counts
            if cds_feature.qualifiers["function"][0] == "unknown function":
                original_functions_count_dict[record_id]["unknown function"] += 1
            elif cds_feature.qualifiers["function"][0] == "transcription regulation":
                original_functions_count_dict[record_id][
                    "transcription regulation"
                ] += 1
            elif cds_feature.qualifiers["function"][0] == "tail":
                original_functions_count_dict[record_id]["tail"] += 1
            elif cds_feature.qualifiers["function"][0] == "other":
                original_functions_count_dict[record_id]["other"] += 1
            elif cds_feature.qualifiers["function"][0] == "moron":
                original_functions_count_dict[record_id][
                    "moron, auxiliary metabolic gene and host takeover"
                ] += 1
            elif cds_feature.qualifiers["function"][0] == "lysis":
                original_functions_count_dict[record_id]["lysis"] += 1
            elif cds_feature.qualifiers["function"][0] == "integration and excision":
                original_functions_count_dict[record_id][
                    "integration and excision"
                ] += 1
            elif cds_feature.qualifiers["function"][0] == "head and packaging":
                original_functions_count_dict[record_id]["head and packaging"] += 1
            elif cds_feature.qualifiers["function"][0] == "DNA":
                original_functions_count_dict[record_id][
                    "DNA, RNA and nucleotide metabolism"
                ] += 1
            elif cds_feature.qualifiers["function"][0] == "connector":
                original_functions_count_dict[record_id]["connector"] += 1

            # now the updated dictionary
            # If record_id does not exist in result_dict, an empty dictionary {} is returned as the default value.
            # prevents KeyError
            if cds_id in result_dict.get(record_id, {}):
                # if cds_id in result_dict[record_id].keys():
                # increase the phrog count
                new_functions_count_dict[record_id]["phrog_count"] += 1
                combined_functions_count_dict[record_id]["phrog_count"] += 1
                # update the counts
                if result_dict[record_id][cds_id]["function"] == "unknown function":
                    new_functions_count_dict[record_id]["unknown function"] += 1
                    combined_functions_count_dict[record_id]["unknown function"] += 1
                elif (
                    result_dict[record_id][cds_id]["function"]
                    == "transcription regulation"
                ):
                    new_functions_count_dict[record_id]["transcription regulation"] += 1
                    combined_functions_count_dict[record_id][
                        "transcription regulation"
                    ] += 1
                elif result_dict[record_id][cds_id]["function"] == "tail":
                    new_functions_count_dict[record_id]["tail"] += 1
                    combined_functions_count_dict[record_id]["tail"] += 1
                elif result_dict[record_id][cds_id]["function"] == "other":
                    new_functions_count_dict[record_id]["other"] += 1
                    combined_functions_count_dict[record_id]["other"] += 1
                elif (
                    result_dict[record_id][cds_id]["function"]
                    == "moron, auxiliary metabolic gene and host takeover"
                ):
                    new_functions_count_dict[record_id][
                        "moron, auxiliary metabolic gene and host takeover"
                    ] += 1
                    combined_functions_count_dict[record_id][
                        "moron, auxiliary metabolic gene and host takeover"
                    ] += 1
                elif result_dict[record_id][cds_id]["function"] == "lysis":
                    new_functions_count_dict[record_id]["lysis"] += 1
                    combined_functions_count_dict[record_id]["lysis"] += 1
                elif (
                    result_dict[record_id][cds_id]["function"]
                    == "integration and excision"
                ):
                    new_functions_count_dict[record_id]["integration and excision"] += 1
                    combined_functions_count_dict[record_id][
                        "integration and excision"
                    ] += 1
                elif result_dict[record_id][cds_id]["function"] == "head and packaging":
                    new_functions_count_dict[record_id]["head and packaging"] += 1
                    combined_functions_count_dict[record_id]["head and packaging"] += 1
                elif (
                    result_dict[record_id][cds_id]["function"]
                    == "DNA, RNA and nucleotide metabolism"
                ):
                    new_functions_count_dict[record_id][
                        "DNA, RNA and nucleotide metabolism"
                    ] += 1
                    combined_functions_count_dict[record_id][
                        "DNA, RNA and nucleotide metabolism"
                    ] += 1
                elif result_dict[record_id][cds_id]["function"] == "connector":
                    new_functions_count_dict[record_id]["connector"] += 1
                    combined_functions_count_dict[record_id]["connector"] += 1

                # update the phrog if different
                # same phrog
                if (
                    result_dict[record_id][cds_id]["phrog"]
                    == cds_feature.qualifiers["phrog"][0]
                ):
                    new_functions_count_dict[record_id]["same_phrogs"] += 1
                # different phrog
                if (
                    result_dict[record_id][cds_id]["phrog"]
                    != cds_feature.qualifiers["phrog"][0]
                ):
                    # where there was no phrog in pharokka
                    if cds_feature.qualifiers["phrog"][0] == "No_PHROG":
                        new_functions_count_dict[record_id]["changed_phrogs"] += 1
                        new_functions_count_dict[record_id]["foldseek_only_phrogs"] += 1
                        updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                            0
                        ] = result_dict[record_id][cds_id]["phrog"]
                        updated_cds_dict[record_id][cds_id].qualifiers["product"][
                            0
                        ] = result_dict[record_id][cds_id]["product"]
                        updated_cds_dict[record_id][cds_id].qualifiers["function"][
                            0
                        ] = result_dict[record_id][cds_id]["function"]
                    # different phrog to pharokka
                    else:
                        # if the foldseek result isn't have unknown function then update
                        if (
                            result_dict[record_id][cds_id]["function"]
                            != "unknown function"
                        ):
                            new_functions_count_dict[record_id]["changed_phrogs"] += 1
                            new_functions_count_dict[record_id]["different_phrogs"] += 1
                            # update
                            updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                0
                            ] = result_dict[record_id][cds_id]["phrog"]
                            updated_cds_dict[record_id][cds_id].qualifiers["product"][
                                0
                            ] = result_dict[record_id][cds_id]["product"]
                            updated_cds_dict[record_id][cds_id].qualifiers["function"][
                                0
                            ] = result_dict[record_id][cds_id]["function"]
            else:
                # will not be in results dict - therefore in foldseek dict will be unknown function function
                new_functions_count_dict[record_id]["unknown function"] += 1
                if cds_feature.qualifiers["phrog"][0] != "No_PHROG":
                    new_functions_count_dict[record_id]["pharokka_only_phrogs"] += 1
                    combined_functions_count_dict[record_id]["phrog_count"] += 1
                    if cds_feature.qualifiers["function"][0] == "unknown function":
                        combined_functions_count_dict[record_id][
                            "unknown function"
                        ] += 1
                    elif (
                        cds_feature.qualifiers["function"][0]
                        == "transcription regulation"
                    ):
                        combined_functions_count_dict[record_id][
                            "transcription regulation"
                        ] += 1
                    elif cds_feature.qualifiers["function"][0] == "tail":
                        combined_functions_count_dict[record_id]["tail"] += 1
                    elif cds_feature.qualifiers["function"][0] == "other":
                        combined_functions_count_dict[record_id]["other"] += 1
                    elif cds_feature.qualifiers["function"][0] == "moron":
                        combined_functions_count_dict[record_id][
                            "moron, auxiliary metabolic gene and host takeover"
                        ] += 1
                    elif cds_feature.qualifiers["function"][0] == "lysis":
                        combined_functions_count_dict[record_id]["lysis"] += 1
                    elif (
                        cds_feature.qualifiers["function"][0]
                        == "integration and excision"
                    ):
                        combined_functions_count_dict[record_id][
                            "integration and excision"
                        ] += 1
                    elif cds_feature.qualifiers["function"][0] == "head and packaging":
                        combined_functions_count_dict[record_id][
                            "head and packaging"
                        ] += 1
                    elif cds_feature.qualifiers["function"][0] == "DNA":
                        combined_functions_count_dict[record_id][
                            "DNA, RNA and nucleotide metabolism"
                        ] += 1
                    elif cds_feature.qualifiers["function"][0] == "connector":
                        combined_functions_count_dict[record_id]["connector"] += 1
                else:
                    # no hits in either
                    combined_functions_count_dict[record_id]["unknown function"] += 1

    # Convert the nested dictionary to a Pandas DataFrame
    pharokka_df = pd.DataFrame.from_dict(original_functions_count_dict, orient="index")
    pharokka_df["contig_id"] = pharokka_df.index
    pharokka_df = pharokka_df[
        ["contig_id"] + [col for col in pharokka_df.columns if col != "contig_id"]
    ]
    pharokka_tsv: Path = Path(output) / "pharokka_functions_output.tsv"
    pharokka_df.to_csv(pharokka_tsv, sep="\t", index=False)

    foldseek_df = pd.DataFrame.from_dict(new_functions_count_dict, orient="index")
    foldseek_df["contig_id"] = foldseek_df.index
    foldseek_df = foldseek_df[
        ["contig_id"] + [col for col in foldseek_df.columns if col != "contig_id"]
    ]
    foldseek_tsv: Path = Path(output) / "foldseek_functions_output.tsv"
    foldseek_df.to_csv(foldseek_tsv, sep="\t", index=False)

    combined_df = pd.DataFrame.from_dict(combined_functions_count_dict, orient="index")
    combined_df["contig_id"] = combined_df.index
    combined_df = combined_df[
        ["contig_id"] + [col for col in combined_df.columns if col != "contig_id"]
    ]
    combined_tsv: Path = Path(output) / "combined_functions_output.tsv"
    combined_df.to_csv(combined_tsv, sep="\t", index=False)
