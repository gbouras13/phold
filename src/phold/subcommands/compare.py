#!/usr/bin/env python3

import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

from phold.features.create_foldseek_db import (
    generate_foldseek_db_from_aa_3di,
    generate_foldseek_db_from_pdbs,
)
from phold.features.run_foldseek import create_result_tsv, run_foldseek_search
from phold.features.split_3Di import split_3di_fasta_by_prob
from phold.io.handle_genbank import write_genbank
from phold.io.sub_db_outputs import create_sub_db_outputs
from phold.results.topfunction import calculate_topfunctions_results, get_topfunctions


def subcommand_compare(
    gb_dict: Dict[str, Dict[str, Union[SeqRecord, SeqFeature]]],
    output: Path,
    threads: int,
    evalue: float,
    card_vfdb_evalue: float,
    sensitivity: float,
    database: Path,
    prefix: str,
    predictions_dir: Optional[Path],
    pdb: bool,
    pdb_dir: Optional[Path],
    logdir: Path,
    filter_pdbs: bool,
    split: bool,
    split_threshold: float,
    remote_flag: bool,
    proteins_flag: bool,
    fasta_flag: bool,
    separate: bool,
    max_seqs: int,
) -> bool:
    """
    Compare 3Di or PDB structures to the Phold DB

    Parameters:
        gb_dict (Dict[str, Dict[str, Union[SeqRecord, SeqFeature]]]): Nested dictionary containing genomic data.
        output (Path): Path to the output directory.
        threads (int): Number of threads to use.
        evalue (float): E-value threshold.
        card_vfdb_evalue (float): E-value threshold for CARD and VFDB databases.
        sensitivity (float): Sensitivity threshold.
        database (Path): Path to the reference database.
        prefix (str): Prefix for output files.
        predictions_dir (Optional[Path]): Path to the directory containing predictions.
        pdb (bool): Flag indicating whether PDB files are used.
        pdb_dir (Optional[Path]): Path to the directory containing PDB files.
        logdir (Path): Path to the directory for log files.
        filter_pdbs (bool): Flag indicating whether to filter PDB files.
        split (bool): Flag indicating whether to split the DBs and run foldseek twice.
        split_threshold (float): Threshold for splitting the DBs  and run foldseek twice.
        remote_flag (bool): Flag indicating whether the analysis is remote.
        proteins_flag (bool): Flag indicating whether proteins are used.
        fasta_flag (bool): Flag indicating whether FASTA files are used.
        separate (bool): Flag indicating whether to separate the analysis.
        max_seqs (int): Maximum results per query sequence allowed to pass the prefilter for foldseek.

    Returns:
        bool: True if sub-databases are created successfully, False otherwise.
    """

    if predictions_dir is None and pdb is False:
        logger.error(f"You did not specify --predictions_dir or --pdb. Please check ")

    if proteins_flag is True:
        cds_dict = gb_dict
        non_cds_dict = {}
    else:
        # Create a nested dictionary to store CDS features by contig ID
        cds_dict = {}

        # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
        for record_id, record in gb_dict.items():
            cds_dict[record_id] = {}

            for cds_feature in record.features:
                # for all pharokka genbanks, correct the bad phrog cats
                if fasta_flag is False:
                    if cds_feature.type == "CDS":
                        # update DNA, RNA and nucleotide metabolism from pharokka as it is broken as of 1.6.1
                        if cds_feature.qualifiers["function"][0] == "DNA":
                            cds_feature.qualifiers["function"][
                                0
                            ] = "DNA, RNA and nucleotide metabolism"
                            cds_feature.qualifiers["function"] = [
                                cds_feature.qualifiers["function"][0]
                            ]  # Keep only the first element
                        # moron, auxiliary metabolic gene and host takeover as it is broken as of 1.6.1
                        if cds_feature.qualifiers["function"][0] == "moron":
                            cds_feature.qualifiers["function"][
                                0
                            ] = "moron, auxiliary metabolic gene and host takeover"
                            cds_feature.qualifiers["function"] = [
                                cds_feature.qualifiers["function"][0]
                            ]  # Keep only the first element

                        # update No_PHROGs_HMM to No_PHROGs - if input is from pharokka --fast
                        if cds_feature.qualifiers["phrog"][0] == "No_PHROGs_HMM":
                            cds_feature.qualifiers["phrog"][0] = "No_PHROG"

                        cds_dict[record_id][
                            cds_feature.qualifiers["ID"][0]
                        ] = cds_feature
                # for fasta, will be fine to just add
                else:
                    cds_dict[record_id][cds_feature.qualifiers["ID"]] = cds_feature

        # non cds dict for later in genbank
        non_cds_dict = {}

        # makes the nested dictionary {contig_id:{non_cds_id: non_cds_feature}}
        for record_id, record in gb_dict.items():
            non_cds_dict[record_id] = {}

            for non_cds_feature in record.features:
                if non_cds_feature.type != "CDS":
                    non_cds_dict[record_id][
                        non_cds_feature.qualifiers["ID"][0]
                    ] = non_cds_feature

    # input predictions or structures
    if pdb is False:
        # prostT5
        fasta_aa_input: Path = Path(predictions_dir) / f"{prefix}_aa.fasta"
        fasta_3di_input: Path = Path(predictions_dir) / f"{prefix}_3di.fasta"

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"

    ## copy the AA and 3Di from predictions directory if pdb is false and phold compare is the command
    if pdb is False:
        # if remote, these will not exist
        if remote_flag is False:
            if fasta_3di_input.exists():
                logger.info(
                    f"Checked that the 3Di CDS file {fasta_3di_input} exists from phold predict"
                )
                shutil.copyfile(fasta_3di_input, fasta_3di)
            else:
                logger.error(
                    f"The 3Di CDS file {fasta_3di_input} does not exist. Please run phold predict and/or check the prediction directory {predictions_dir}"
                )
            # copy the aa to file
            if fasta_aa_input.exists():
                logger.info(
                    f"Checked that the AA CDS file {fasta_aa_input} exists from phold predict."
                )
                shutil.copyfile(fasta_aa_input, fasta_aa)
            else:
                logger.error(
                    f"The AA CDS file {fasta_aa_input} does not exist. Please run phold predict and/or check the prediction directory {predictions_dir}"
                )
    ## write the AAs to file if pdb is true
    else:
        ## write the CDS to file
        logger.info(f"Writing the AAs to file {fasta_aa}.")
        with open(fasta_aa, "w+") as out_f:
            for record_id, rest in cds_dict.items():
                aa_contig_dict = cds_dict[record_id]

                # writes the CDS to file
                for seq_id, cds_feature in aa_contig_dict.items():
                    out_f.write(f">{record_id}:{seq_id}\n")
                    out_f.write(f"{cds_feature.qualifiers['translation'][0]}\n")

    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    if pdb is True:
        logger.info("Creating a foldseek query database from structures.")

        filtered_pdbs_path: Path = Path(output) / "filtered_pdbs"
        if filter_pdbs is True:
            logger.info(
                f"--filter_pdbs is {filter_pdbs}. Therefore only the .pdb structure files with matching CDS ids will be copied and compared."
            )
            filtered_pdbs_path.mkdir(parents=True, exist_ok=True)

        generate_foldseek_db_from_pdbs(
            fasta_aa,
            foldseek_query_db_path,
            pdb_dir,
            filtered_pdbs_path,
            logdir,
            prefix,
            filter_pdbs,
        )
    else:
        # split
        if split is True:
            probs_3di: Path = (
                Path(predictions_dir) / f"{prefix}_prostT5_3di_mean_probabilities.csv"
            )
            split_3di_fasta_by_prob(
                fasta_aa, fasta_3di, probs_3di, output, split_threshold
            )
            logger.info(
                f"--split is {split} with threshold {split_threshold}. Generating Foldseek query DBs for high and low ProstT5 probability subsets."
            )

            high_prob_fasta_aa_out_path: Path = Path(output) / "high_prob_aa.fasta"
            low_prob_fasta_aa_out_path: Path = Path(output) / "low_prob_aa.fasta"
            high_prob_fasta_3di_out_path: Path = Path(output) / "high_prob_3di.fasta"
            low_prob_fasta_3di_out_path: Path = Path(output) / "low_prob_3di.fasta"

            # high
            if "high_prostt5_prob" == prefix:
                logger.error(f"Please choose a different {prefix}. ")

            generate_foldseek_db_from_aa_3di(
                high_prob_fasta_aa_out_path,
                high_prob_fasta_3di_out_path,
                foldseek_query_db_path,
                logdir,
                "high_prostt5_prob",
            )

            # low
            if "low_prostt5_prob" == prefix:
                logger.error(f"Please choose a different {prefix}.t ")
            generate_foldseek_db_from_aa_3di(
                low_prob_fasta_aa_out_path,
                low_prob_fasta_3di_out_path,
                foldseek_query_db_path,
                logdir,
                "low_prostt5_prob",
            )
        # default (just prostt5)
        else:
            generate_foldseek_db_from_aa_3di(
                fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
            )

    #####
    # foldseek search
    #####

    if split is True:

        #############
        # high - vs structures
        #############
        query_db: Path = Path(foldseek_query_db_path) / "high_prostt5_prob"
        target_db: Path = Path(database) / "all_phold_structures"

        # make result and temp dirs
        result_db_base: Path = Path(output) / "result_db"
        result_db_base.mkdir(parents=True, exist_ok=True)
        result_db: Path = Path(result_db_base) / "result_db_high"

        temp_db: Path = Path(output) / "temp_db"
        temp_db.mkdir(parents=True, exist_ok=True)

        # run foldseek search
        run_foldseek_search(
            query_db,
            target_db,
            result_db,
            temp_db,
            threads,
            logdir,
            evalue,
            sensitivity,
            max_seqs,
        )

        # make result tsv for high prob vs structure db
        result_high_tsv: Path = Path(output) / "foldseek_results_high.tsv"
        create_result_tsv(query_db, target_db, result_db, result_high_tsv, logdir)

        #############
        # low - vs prostt5 db
        ############

        query_db: Path = Path(foldseek_query_db_path) / "low_prostt5_prob"
        target_db: Path = Path(database) / "all_phold_prostt5"

        # make result and temp dirs
        result_db_base: Path = Path(output) / "result_db"
        result_db_base.mkdir(parents=True, exist_ok=True)
        result_db: Path = Path(result_db_base) / "result_db_low"

        temp_db: Path = Path(output) / "temp_db"
        temp_db.mkdir(parents=True, exist_ok=True)

        # run foldseek search
        run_foldseek_search(
            query_db,
            target_db,
            result_db,
            temp_db,
            threads,
            logdir,
            evalue,
            sensitivity,
            max_seqs,
        )

        # make result tsv
        result_low_tsv: Path = Path(output) / "foldseek_results_low.tsv"
        create_result_tsv(query_db, target_db, result_db, result_low_tsv, logdir)

        # create combined tsv
        high_df = pd.read_csv(result_high_tsv, sep="\t", header=None)
        low_df = pd.read_csv(result_low_tsv, sep="\t", header=None)

        # combined
        combined_df = pd.concat([high_df, low_df])
        result_tsv: Path = Path(output) / "foldseek_results.tsv"
        combined_df.to_csv(result_tsv, sep="\t", header=False, index=False)

    # no split
    else:

        short_db_name = prefix
        database_name = "all_phold_structures"
        if short_db_name == database_name:
            logger.error(
                f"Please choose a different {prefix} as this conflicts with the {database_name}"
            )

        query_db: Path = Path(foldseek_query_db_path) / short_db_name
        target_db: Path = Path(database) / database_name

        # make result and temp dirs
        result_db_base: Path = Path(output) / "result_db"
        result_db_base.mkdir(parents=True, exist_ok=True)
        result_db: Path = Path(result_db_base) / "result_db"

        temp_db: Path = Path(output) / "temp_db"
        temp_db.mkdir(parents=True, exist_ok=True)

        # run foldseek search
        run_foldseek_search(
            query_db,
            target_db,
            result_db,
            temp_db,
            threads,
            logdir,
            evalue,
            sensitivity,
            max_seqs,
        )

        # make result tsv
        result_tsv: Path = Path(output) / "foldseek_results.tsv"
        create_result_tsv(query_db, target_db, result_db, result_tsv, logdir)

    ########
    # get topfunction
    ########

    result_tsv: Path = Path(output) / "foldseek_results.tsv"

    # get top hit with non-unknown function for each CDS
    # also calculate the weighted bitscore df

    filtered_topfunctions_df, weighted_bitscore_df = get_topfunctions(
        result_tsv, database, database_name, pdb, card_vfdb_evalue, proteins_flag
    )

    # update the CDS dictionary with the tophits
    updated_cds_dict, filtered_tophits_df = calculate_topfunctions_results(
        filtered_topfunctions_df, cds_dict, output, pdb, proteins_flag, fasta_flag
    )

    # generate per CDS foldseek information df and write to genbank
    per_cds_df = write_genbank(
        updated_cds_dict,
        non_cds_dict,
        prefix,
        gb_dict,
        output,
        proteins_flag,
        separate,
        fasta_flag,
    )

    # if prostt5, query will repeat contig_id in query - convert to cds_id
    if pdb is False and proteins_flag is False:
        weighted_bitscore_df[["contig_id", "cds_id"]] = weighted_bitscore_df[
            "query"
        ].str.split(":", expand=True, n=1)
        weighted_bitscore_df = weighted_bitscore_df.drop(columns=["query", "contig_id"])

    # otherwise for pdb or proteins query will just be the cds_id so rename
    else:
        weighted_bitscore_df.rename(columns={"query": "cds_id"}, inplace=True)

    if proteins_flag is False:
        columns_to_drop = ["contig_id", "phrog", "product", "function"]
    else:
        columns_to_drop = ["phrog", "product", "function"]
    # drop dupe columns
    filtered_tophits_df = filtered_tophits_df.drop(columns=columns_to_drop)

    merged_df = per_cds_df.merge(filtered_tophits_df, on="cds_id", how="left")
    merged_df = merged_df.merge(weighted_bitscore_df, on="cds_id", how="left")

    ########
    #### add annotation source
    ########

    # Define a function to apply to each row to determine the annotation source
    def determine_annotation_source(row):
        if row["phrog"] == "No_PHROG":
            return "none"
        elif pd.isnull(row["bitscore"]):
            return "pharokka"
        else:
            return "foldseek"

    # Create a new column order with 'annotation_method' moved after 'product'
    merged_df["annotation_method"] = merged_df.apply(
        determine_annotation_source, axis=1
    )

    product_index = merged_df.columns.get_loc("product")

    new_column_order = (
        list(merged_df.columns[: product_index + 1])
        + ["annotation_method"]
        + list(merged_df.columns[product_index + 1 : -1])
    )
    merged_df = merged_df.reindex(columns=new_column_order)

    merged_df_path: Path = Path(output) / f"{prefix}_per_cds_predictions.tsv"
    merged_df.to_csv(merged_df_path, index=False, sep="\t")

    # save vfdb card acr defensefinder hits with more metadata
    sub_dbs_created = create_sub_db_outputs(merged_df, database, output)

    # save the function counts is not proteins

    if proteins_flag is False:

        contig_ids = merged_df["contig_id"].unique()

        # get list of all functions counts
        functions_list = []

        for contig in contig_ids:

            contig_df = merged_df[merged_df["contig_id"] == contig]

            cds_count = len(contig_df)
            # get counts of functions and cds
            # all 10 PHROGs categories
            connector_count = len(contig_df[contig_df["function"] == "connector"])
            metabolism_count = len(
                contig_df[contig_df["function"] == "DNA, RNA and nucleotide metabolism"]
            )
            head_count = len(contig_df[contig_df["function"] == "head and packaging"])
            integration_count = len(
                contig_df[contig_df["function"] == "integration and excision"]
            )
            lysis_count = len(contig_df[contig_df["function"] == "lysis"])
            moron_count = len(
                contig_df[
                    contig_df["function"]
                    == "moron, auxiliary metabolic gene and host takeover"
                ]
            )
            other_count = len(contig_df[contig_df["function"] == "other"])
            tail_count = len(contig_df[contig_df["function"] == "tail"])
            transcription_count = len(
                contig_df[contig_df["function"] == "transcription regulation"]
            )
            unknown_count = len(contig_df[contig_df["function"] == "unknown function"])

            acr_count = len(contig_df[contig_df["phrog"] == "acr"])

            vfdb_count = len(contig_df[contig_df["phrog"] == "vfdb"])

            card_count = len(contig_df[contig_df["phrog"] == "card"])

            defensefinder_count = len(contig_df[contig_df["phrog"] == "defensefinder"])

            # create count list  for the dataframe
            count_list = [
                cds_count,
                connector_count,
                metabolism_count,
                head_count,
                integration_count,
                lysis_count,
                moron_count,
                other_count,
                tail_count,
                transcription_count,
                unknown_count,
            ]

            description_list = [
                "CDS",
                "connector",
                "DNA, RNA and nucleotide metabolism",
                "head and packaging",
                "integration and excision",
                "lysis",
                "moron, auxiliary metabolic gene and host takeover",
                "other",
                "tail",
                "transcription regulation",
                "unknown function",
            ]
            contig_list = [
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
            ]
            # cds df
            cds_df = pd.DataFrame(
                {
                    "Description": description_list,
                    "Count": count_list,
                    "Contig": contig_list,
                }
            )

            vfdb_row = pd.DataFrame(
                {
                    "Description": ["VFDB_Virulence_Factors"],
                    "Count": [vfdb_count],
                    "Contig": [contig],
                }
            )
            card_row = pd.DataFrame(
                {
                    "Description": ["CARD_AMR"],
                    "Count": [card_count],
                    "Contig": [contig],
                }
            )

            acr_row = pd.DataFrame(
                {
                    "Description": ["ACR_anti_crispr"],
                    "Count": [acr_count],
                    "Contig": [contig],
                }
            )

            defensefinder_row = pd.DataFrame(
                {
                    "Description": ["Defensefinder"],
                    "Count": [defensefinder_count],
                    "Contig": [contig],
                }
            )

            # eappend it all to functions_list
            functions_list.append(cds_df)
            functions_list.append(vfdb_row)
            functions_list.append(card_row)
            functions_list.append(acr_row)
            functions_list.append(defensefinder_row)

        # combine all contigs into one final df
        description_total_df = pd.concat(functions_list)

        descriptions_total_path: Path = Path(output) / f"{prefix}_all_cds_functions.tsv"
        description_total_df.to_csv(descriptions_total_path, index=False, sep="\t")

    return sub_dbs_created
