#!/usr/bin/env python3
import copy
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from loguru import logger


def get_topfunctions(
    result_tsv: Path, database: Path, database_name: str, pdb: bool
) -> pd.DataFrame:
    logger.info("Processing Foldseek output.")
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

    foldseek_df = pd.read_csv(
        result_tsv, delimiter="\t", index_col=False, names=col_list
    )

    # gets the cds
    if pdb is False:
        # prostt5
        foldseek_df[["contig_id", "cds_id"]] = foldseek_df["query"].str.split(
            ":", expand=True, n=1
        )
    else:
        foldseek_df["cds_id"] = foldseek_df["query"].str.replace(".pdb", "")

    # clean up later
    if (
        database_name == "all_phrogs"
        or database_name == "all_envhogs"
        or database_name == "all_phrogs_pdb"
        or database_name == "all_phold_structures"
        or database_name == "all_phold_prostt5"
    ):
        foldseek_df["target"] = foldseek_df["target"].str.replace(".pdb", "")
        # split the target column as this will have phrog:protein
        foldseek_df[["phrog", "tophit_protein"]] = foldseek_df["target"].str.split(
            ":", expand=True, n=1
        )

    foldseek_df = foldseek_df.drop(columns=["target"])
    foldseek_df["phrog"] = foldseek_df["phrog"].str.replace("phrog_", "")

    mask = foldseek_df["phrog"].str.startswith("envhog_")
    # strip off envhog
    foldseek_df.loc[mask, "phrog"] = foldseek_df.loc[mask, "phrog"].str.replace(
        "envhog_", ""
    )
    # add envhog to protein
    foldseek_df.loc[mask, "tophit_protein"] = (
        "envhog_" + foldseek_df.loc[mask, "tophit_protein"]
    )

    foldseek_df["phrog"] = foldseek_df["phrog"].astype("str")

    # read in the mapping tsv
    phrog_annot_mapping_tsv: Path = Path(database) / "phold_annots.tsv"
    phrog_mapping_df = pd.read_csv(phrog_annot_mapping_tsv, sep="\t")
    phrog_mapping_df["phrog"] = phrog_mapping_df["phrog"].astype("str")

    # join the dfs
    foldseek_df = foldseek_df.merge(phrog_mapping_df, on="phrog", how="left")
    # Replace NaN values in the 'product' column with 'hypothetical protein'
    foldseek_df["product"] = foldseek_df["product"].fillna("hypothetical protein")
    foldseek_df["evalue"] = foldseek_df["evalue"].apply(
        lambda x: "{:.3e}".format(float(x))
    )

    def custom_nsmallest(group):
        # where all the
        if all(group["product"] == "hypothetical protein"):
            min_row_index = group["evalue"].idxmin()
            # Get the entire row
            return group.loc[min_row_index]
        else:
            group = group[group["product"] != "hypothetical protein"]
            min_row_index = group["evalue"].idxmin()
            # Get the entire row
            return group.loc[min_row_index]

    topfunction_df = (
        foldseek_df.groupby("query", group_keys=True)
        .apply(custom_nsmallest)
        .reset_index(drop=True)
    )

    topfunction_dict = dict(zip(topfunction_df["query"], topfunction_df["function"]))

    # Remove the original 'query' column
    topfunction_df = topfunction_df.drop(columns=["query"])

    # scientific notation to 3dp
    topfunction_df["evalue"] = topfunction_df["evalue"].apply(
        lambda x: "{:.3e}".format(float(x))
    )

    def weighted_function(group):

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

        query = group["query"].iloc[0]
        tophit_function = topfunction_dict[query]

        # function_counts = group['function'].value_counts().to_dict()

        # normalise counts by total bitscore
        weighted_counts_normalised = {}
        # total_bitscore = group['bitscore'].sum()
        bitscore_by_function = group.groupby("function")["bitscore"].sum().to_dict()


        total_functional_bitscore = group[group["function"] != "unknown function"]["bitscore"].sum()

        if total_functional_bitscore == 0:
        # if tophit_function == "unknown function":
            top_bitscore_function = "unknown function"
            top_bitscore_perc = 0
        
        # everything except unknown function
        # get total bitscore of the hits with function
        else:
            # total_function_bitscore = group[group["function"] != "unknown function"]["bitscore"].sum()

            # get the weighted bitscore
            for key, value in bitscore_by_function.items():
                if key != "unknown function":
                    weighted_counts_normalised[key] = round(
                        value / total_functional_bitscore, 2
                    )

            top_bitscore_function = max(
                weighted_counts_normalised, key=weighted_counts_normalised.get
            )
            top_bitscore_perc = max(weighted_counts_normalised.values())

        d = {
            "tophit_function": [tophit_function],
            "function_with_highest_bitscore_percentage": [top_bitscore_function],
            "top_bitscore_percentage_not_unknown": [top_bitscore_perc],
            "head_and_packaging_bitscore_percentage": [
                weighted_counts_normalised.get("head and packaging", 0)
            ],
            "integration_and_excision_bitscore_percentage": [
                weighted_counts_normalised.get("integration and excision", 0)
            ],
            "tail bitscore_percentage": [weighted_counts_normalised.get("tail", 0)],
            "moron_auxiliary_metabolic_gene_and_host_takeover_bitscore_percentage": [
                weighted_counts_normalised.get(
                    "moron, auxiliary metabolic gene and host takeover", 0
                )
            ],
            "DNA_RNA_and_nucleotide_metabolism bitscore_percentage": [
                weighted_counts_normalised.get("DNA, RNA and nucleotide metabolism", 0)
            ],
            "connector_bitscore_percentage": [
                weighted_counts_normalised.get("connector", 0)
            ],
            "transcription_regulation_bitscore_percentage": [
                weighted_counts_normalised.get("transcription regulation", 0)
            ],
            "lysis_bitscore_percentage": [weighted_counts_normalised.get("lysis", 0)],
            "other_bitscore_percentage": [weighted_counts_normalised.get("other", 0)],
            "unknown_function_bitscore_percentage": [
                weighted_counts_normalised.get("unknown function", 0)
            ],
        }

        weighted_bitscore_df = pd.DataFrame(data=d)

        return weighted_bitscore_df

    weighted_bitscore_df = foldseek_df.groupby("query", group_keys=True).apply(
        weighted_function
    )

    weighted_bitscore_df.reset_index(inplace=True)
    weighted_bitscore_df["query"] = weighted_bitscore_df["query"].str.replace(
        ".pdb", ""
    )
    weighted_bitscore_df = weighted_bitscore_df.drop(columns=["level_1"])

    return topfunction_df, weighted_bitscore_df


def calculate_topfunctions_results(
    filtered_tophits_df: pd.DataFrame, cds_dict: dict, output: Path, pdb: bool
) -> None:
    # Convert the DataFrame to a nested dictionary
    result_dict = {}

    # so I can match with the df row below
    cds_record_dict = {}

    for record_id, cds_entries in cds_dict.items():
        result_dict[record_id] = {}
        for cds_id, cds_info in cds_entries.items():
            result_dict[record_id][cds_id] = {}
            cds_record_dict[cds_id] = record_id

    # Get record_id for every cds_id and merge
    if pdb is True:
        cds_record_df = pd.DataFrame(
            list(cds_record_dict.items()), columns=["cds_id", "contig_id"]
        )
        filtered_tophits_df = pd.merge(
            filtered_tophits_df, cds_record_df, on="cds_id", how="left"
        )

    # Move "contig_id" and 'cds_id' to the front
    column_order = ["contig_id", "cds_id"] + [
        col
        for col in filtered_tophits_df.columns
        if col != "contig_id" and col != "cds_id"
    ]
    filtered_tophits_df = filtered_tophits_df[column_order]

    # combined_tsv_path: Path = Path(output) / "topfunctions.tsv"
    # filtered_tophits_df.to_csv(combined_tsv_path, sep="\t", index=False)

    # loop over all the foldseek tophits
    for _, row in filtered_tophits_df.iterrows():
        record_id = row["contig_id"]
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
        result_dict[record_id][cds_id] = values_dict

    # get counds
    # copy initial cds_dict
    updated_cds_dict = copy.deepcopy(cds_dict)

    original_functions_count_dict = {}
    new_functions_count_dict = {}
    combined_functions_count_dict = {}
    comparison_count_dict = {}

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

    # iterates over the records
    for record_id, record in updated_cds_dict.items():
        cds_count = len(updated_cds_dict[record_id])

        # instantiate the functions dicts
        original_functions_count_dict[record_id] = {}
        new_functions_count_dict[record_id] = {}
        combined_functions_count_dict[record_id] = {}

        # this one holds same_phrogs, changed_phrogs, different_phrogs,  foldseek_only_phrogs, pharokka_only_phrogs
        comparison_count_dict[record_id] = {}
        comparison_count_dict[record_id]["cds_count"] = cds_count

        # instantiate every CDS dictionary

        original_functions_count_dict = initialize_function_counts_dict(
            record_id, original_functions_count_dict, cds_count
        )
        new_functions_count_dict = initialize_function_counts_dict(
            record_id, new_functions_count_dict, cds_count
        )
        combined_functions_count_dict = initialize_function_counts_dict(
            record_id, combined_functions_count_dict, cds_count
        )

        # the extra ones instnatiate manually
        comparison_count_dict[record_id]["phrog_count"] = 0
        comparison_count_dict[record_id]["same_phrogs"] = 0
        comparison_count_dict[record_id]["changed_phrogs"] = 0
        comparison_count_dict[record_id]["kept_phrogs"] = 0
        comparison_count_dict[record_id]["foldseek_only_phrogs"] = 0
        comparison_count_dict[record_id]["different_phrogs"] = 0
        comparison_count_dict[record_id]["pharokka_only_phrogs"] = 0
        comparison_count_dict[record_id]["neither_phrogs"] = 0
        # function
        comparison_count_dict[record_id]["same_function"] = 0
        comparison_count_dict[record_id]["changed_function"] = 0
        comparison_count_dict[record_id]["foldseek_only_non_hyp_function"] = 0
        comparison_count_dict[record_id]["pharokka_only_non_hyp_function"] = 0

        # iterates over the features
        for cds_id, cds_feature in updated_cds_dict[record_id].items():
            # if pharokka has a phrog
            if cds_feature.qualifiers["phrog"][0] != "No_PHROG":
                original_functions_count_dict[record_id]["phrog_count"] += 1

            # original pharokka phrog categories
            pharokka_phrog_function_category = phrog_function_mapping.get(
                cds_feature.qualifiers["function"][0], None
            )

            if pharokka_phrog_function_category:
                original_functions_count_dict[record_id][
                    pharokka_phrog_function_category
                ] += 1

            # if the result_dict is not empty
            if result_dict[record_id][cds_id] != {}:
                # update the combined and foldseek functions counts
                new_functions_count_dict[record_id]["phrog_count"] += 1
                combined_functions_count_dict[record_id]["phrog_count"] += 1
                comparison_count_dict[record_id]["phrog_count"] += 1

                # get the foldseek function
                for function_key in phrog_function_mapping.values():
                    new_functions_count_dict[record_id].setdefault(function_key, 0)
                    combined_functions_count_dict[record_id].setdefault(function_key, 0)

                # function will be None if there is no foldseek hit
                function = result_dict[record_id][cds_id].get("function", None)
                if function is not None:
                    function_key = phrog_function_mapping.get(
                        function, "unknown function"
                    )
                else:
                    function_key = "unknown function"

                # update with new function
                new_functions_count_dict[record_id][function_key] += 1
                combined_functions_count_dict[record_id][function_key] += 1

                # update the phrog if different

                foldseek_phrog = result_dict[record_id][cds_id].get("phrog", None)

                # same phrog as pharokka
                if foldseek_phrog == cds_feature.qualifiers["phrog"][0]:
                    comparison_count_dict[record_id]["same_phrogs"] += 1
                    comparison_count_dict[record_id]["same_function"] += 1
                # different phrog as pharokka
                else:
                    # where there was no phrog in pharokka
                    if cds_feature.qualifiers["phrog"][0] == "No_PHROG":
                        comparison_count_dict[record_id]["changed_phrogs"] += 1
                        comparison_count_dict[record_id]["foldseek_only_phrogs"] += 1
                        updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                            0
                        ] = result_dict[record_id][cds_id]["phrog"]
                        updated_cds_dict[record_id][cds_id].qualifiers["product"][
                            0
                        ] = result_dict[record_id][cds_id]["product"]
                        updated_cds_dict[record_id][cds_id].qualifiers["function"][
                            0
                        ] = result_dict[record_id][cds_id]["function"]

                        # whether function changed
                        if (
                            result_dict[record_id][cds_id]["function"]
                            != "unknown function"
                        ):
                            comparison_count_dict[record_id][
                                "foldseek_only_non_hyp_function"
                            ] += 1
                            comparison_count_dict[record_id]["changed_function"] += 1
                        # otherwise same function
                        else:
                            comparison_count_dict[record_id]["same_function"] += 1

                    # different phrog to pharokka
                    # maybe look at this differently?
                    else:
                        # if the foldseek result doesn't have unknown function then update
                        if (
                            result_dict[record_id][cds_id]["function"]
                            != "unknown function"
                        ):
                            # pharokka hyp, foldseek not
                            if (
                                cds_feature.qualifiers["function"][0]
                                == "unknown function"
                                and result_dict[record_id][cds_id]["function"]
                                != "unknown function"
                            ):
                                comparison_count_dict[record_id][
                                    "foldseek_only_non_hyp_function"
                                ] += 1

                            if (
                                cds_feature.qualifiers["function"][0]
                                != result_dict[record_id][cds_id]["function"]
                            ):
                                comparison_count_dict[record_id][
                                    "changed_function"
                                ] += 1
                            else:
                                comparison_count_dict[record_id]["same_function"] += 1

                            comparison_count_dict[record_id]["changed_phrogs"] += 1
                            comparison_count_dict[record_id]["different_phrogs"] += 1
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
                        # foldseek result has unknown function
                        else:
                            # foldseek result is unknown function and pharokka has known function
                            # keep the pharokka phrogs
                            if (
                                cds_feature.qualifiers["function"][0]
                                != "unknown function"
                            ):
                                comparison_count_dict[record_id]["kept_phrogs"] += 1
                                comparison_count_dict[record_id][
                                    "pharokka_only_non_hyp_function"
                                ] += 1
                            # foldseek result is unknown function and pharokka unknown function
                            # update to foldseek
                            # arguably can change this - will need  comparison_count_dict[record_id]["kept_phrogs"] += 1 if so
                            else:
                                comparison_count_dict[record_id]["changed_phrogs"] += 1
                                updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                    0
                                ] = result_dict[record_id][cds_id]["phrog"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "product"
                                ][0] = result_dict[record_id][cds_id]["product"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "function"
                                ][0] = result_dict[record_id][cds_id]["function"]

                                comparison_count_dict[record_id]["same_function"] += 1

            # no foldseek hits - empty dict
            # will not be in results dict - therefore in foldseek dict will be function function
            else:
                new_functions_count_dict[record_id]["unknown function"] += 1
                # if has a pharokka hit
                if cds_feature.qualifiers["phrog"][0] != "No_PHROG":
                    combined_functions_count_dict[record_id]["phrog_count"] += 1
                    comparison_count_dict[record_id]["phrog_count"] += 1
                    comparison_count_dict[record_id]["pharokka_only_phrogs"] += 1
                    # comparison_count_dict[record_id]["changed_phrogs"] += 1

                    function_value = cds_feature.qualifiers["function"][0]
                    function_key = phrog_function_mapping.get(
                        function_value, "unknown function"
                    )

                    combined_functions_count_dict[record_id].setdefault(function_key, 0)
                    combined_functions_count_dict[record_id][function_key] += 1

                    if (
                        updated_cds_dict[record_id][cds_id].qualifiers["function"][0]
                        != "unknown function"
                    ):
                        comparison_count_dict[record_id][
                            "pharokka_only_non_hyp_function"
                        ] += 1
                    else:
                        comparison_count_dict[record_id]["same_function"] += 1

                # no hits in either
                else:
                    combined_functions_count_dict[record_id]["unknown function"] += 1
                    comparison_count_dict[record_id]["neither_phrogs"] += 1
                    comparison_count_dict[record_id]["same_function"] += 1

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

    comparison_df = pd.DataFrame.from_dict(comparison_count_dict, orient="index")
    comparison_df["contig_id"] = comparison_df.index
    comparison_df = comparison_df[
        ["contig_id"] + [col for col in comparison_df.columns if col != "contig_id"]
    ]
    comparison_tsv: Path = Path(output) / "comparison_output.tsv"
    comparison_df.to_csv(comparison_tsv, sep="\t", index=False)

    return updated_cds_dict, filtered_tophits_df


def initialize_function_counts_dict(record_id, count_dict, cds_count):
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
