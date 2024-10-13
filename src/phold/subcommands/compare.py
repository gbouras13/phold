#!/usr/bin/env python3

import shutil
from pathlib import Path
from typing import Dict, Optional, Union

import polars as pl
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

from phold.features.create_foldseek_db import (
    generate_foldseek_db_from_aa_3di, generate_foldseek_db_from_structures)
from phold.features.run_foldseek import create_result_tsv, run_foldseek_search
from phold.io.handle_genbank import write_genbank
from phold.io.sub_db_outputs import create_sub_db_outputs
from phold.results.topfunction import (calculate_topfunctions_results,
                                       get_topfunctions)


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
    structures: bool,
    structure_dir: Optional[Path],
    logdir: Path,
    filter_structures: bool,
    remote_flag: bool,
    proteins_flag: bool,
    fasta_flag: bool,
    separate: bool,
    max_seqs: int,
    only_representatives: bool,
    ultra_sensitive: bool
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
        structures (bool): Flag indicating whether structures files are used.
        structure_dir (Optional[Path]): Path to the directory containing structures (.pdb or .cif) files.
        logdir (Path): Path to the directory for log files.
        filter_structuress (bool): Flag indicating whether to filter structure files.
        remote_flag (bool): Flag indicating whether the analysis is remote.
        proteins_flag (bool): Flag indicating whether proteins are used.
        fasta_flag (bool): Flag indicating whether FASTA files are used.
        separate (bool): Flag indicating whether to separate the analysis.
        max_seqs (int): Maximum results per query sequence allowed to pass the prefilter for foldseek.
        only_representatives (bool): Whether to search against representatives only (turn off --cluster-search 1)
        ultra_sensitive (bool): Whether to skip foldseek prefilter for maximum sensitivity

    Returns:
        bool: True if sub-databases are created successfully, False otherwise.
    """

    if predictions_dir is None and structures is False:
        logger.error(
            f"You did not specify --predictions_dir or --structures. Please check "
        )

    if structures and structure_dir is None:
        logger.error(
            f"You specified --structures but you did not specify --structure_dir. Please check "
        )

    if structure_dir and structures is False:
        logger.error(
            f"You specified --structure_dir but you did not specify --structures. Please check "
        )

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
                        # first try Pharokka (ID) - do the updating
                        try:
                            # cds_id = cds_feature.qualifiers["ID"][0]
                            # update DNA, RNA and nucleotide metabolism from pharokka as it is broken as of 1.6.1
                            if "DNA" in cds_feature.qualifiers["function"][0]:
                                cds_feature.qualifiers["function"][
                                    0
                                ] = "DNA, RNA and nucleotide metabolism"
                                cds_feature.qualifiers["function"] = [
                                    cds_feature.qualifiers["function"][0]
                                ]  # Keep only the first element
                            # moron, auxiliary metabolic gene and host takeover as it is broken as of 1.6.1
                            if "moron" in cds_feature.qualifiers["function"][0]:
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

                        # not pharokka - must be from genbank (supported only)
                        except:
                            try:
                                # add these extra fields to make it all play nice
                                cds_feature.qualifiers["ID"] = cds_feature.qualifiers[
                                    "protein_id"
                                ]
                                cds_feature.qualifiers["function"] = []
                                cds_feature.qualifiers["function"].append(
                                    "unknown function"
                                )
                                cds_feature.qualifiers["phrog"] = []
                                cds_feature.qualifiers["phrog"].append("No_PHROG")

                                cds_dict[record_id][
                                    cds_feature.qualifiers["ID"][0]
                                ] = cds_feature

                            except:
                                logger.error(
                                    f"Feature {cds_feature} has no 'ID' or 'protein_id' qualifier in the Genbank file. Please add one in."
                                )

                # for fasta, will be fine to just add
                else:
                    cds_dict[record_id][cds_feature.qualifiers["ID"]] = cds_feature

        # non cds dict for later in genbank
        non_cds_dict = {}

        # makes the nested dictionary {contig_id:{non_cds_id: non_cds_feature}}
        for record_id, record in gb_dict.items():
            non_cds_dict[record_id] = {}

            i = 1
            for non_cds_feature in record.features:
                if non_cds_feature.type != "CDS":
                    try:
                        non_cds_dict[record_id][
                            non_cds_feature.qualifiers["ID"][0]
                        ] = non_cds_feature
                    except:
                        non_cds_dict[record_id][
                            f"non_cds_feature_{i}"
                        ] = non_cds_feature
                        i += 1

    # input predictions or structures
    if structures is False:
        # prostT5
        fasta_aa_input: Path = Path(predictions_dir) / f"{prefix}_aa.fasta"
        fasta_3di_input: Path = Path(predictions_dir) / f"{prefix}_3di.fasta"

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
    fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"

    ## copy the AA and 3Di from predictions directory if structures is false and phold compare is the command
    if structures is False:
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
    ## write the AAs to file if structures is true because can't just copy from prediction_dir
    else:
        ## write the CDS to file
        logger.info(f"Writing the AAs to file {fasta_aa}.")
        with open(fasta_aa, "w+") as out_f:
            for record_id, rest in cds_dict.items():
                aa_contig_dict = cds_dict[record_id]

                # writes the CDS to file
                for seq_id, cds_feature in aa_contig_dict.items():
                    # if proteins, don't want the 'proteins:' as CDS id
                    if proteins_flag:
                        header = f">{seq_id}\n"
                        seq = f"{cds_feature.qualifiers['translation']}\n"
                    else:  # if genbank entry need to take the first seq as it is parsed as a list
                        header = f">{record_id}:{seq_id}\n"
                        seq = f"{cds_feature.qualifiers['translation'][0]}\n"
                    out_f.write(header)
                    out_f.write(seq)

    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    if structures is True:
        logger.info("Creating a foldseek query database from structures.")

        filtered_structures_path: Path = Path(output) / "filtered_structures"
        if filter_structures is True:
            logger.info(
                f"--filter_structures is {filter_structures}. Therefore only the .pdb or .cif structure files with matching CDS ids will be copied and compared."
            )
            filtered_structures_path.mkdir(parents=True, exist_ok=True)

        generate_foldseek_db_from_structures(
            fasta_aa,
            foldseek_query_db_path,
            structure_dir,
            filtered_structures_path,
            logdir,
            prefix,
            filter_structures,
            proteins_flag,
        )
    else:
        generate_foldseek_db_from_aa_3di(
            fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
        )

    short_db_name = prefix

    # # clustered db search
    # if cluster_db is True:
    #     database_name = "all_phold_structures_clustered_searchDB"
    # else:
    #     database_name = "all_phold_structures"

    # clustered DB
    database_name = "all_phold_structures_clustered_searchDB"

    if short_db_name == database_name:
        logger.error(
            f"Please choose a different {prefix} as this conflicts with the {database_name}"
        )

    #####
    # foldseek search
    #####

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
        only_representatives,
        ultra_sensitive,
    )

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results.tsv"

    # target_db is all_phold_structures regardless of the clustered search mode
    # needs all_phold_structures and all_phold_structures_h
    # delete the rest to save some space

    target_db: Path = Path(database) / "all_phold_structures"
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir)

    ########
    # get topfunction
    ########

    result_tsv: Path = Path(output) / "foldseek_results.tsv"

    # Get top hit with non-unknown function for each CDS and calculate the weighted bitscore
    filtered_topfunctions_df, weighted_bitscore_df = get_topfunctions(
        result_tsv, database, database_name, structures, card_vfdb_evalue, proteins_flag
    )

    # Update the CDS dictionary with the top hits
    updated_cds_dict, filtered_tophits_df, source_dict = calculate_topfunctions_results(
        filtered_topfunctions_df,
        cds_dict,
        output,
        structures,
        proteins_flag,
        fasta_flag,
    )

    # Generate per CDS foldseek information and write to GenBank
    per_cds_df = write_genbank(
        updated_cds_dict,
        non_cds_dict,
        source_dict,
        prefix,
        gb_dict,
        output,
        proteins_flag,
        separate,
        fasta_flag,
    )

    # Handle query modifications based on conditions
    if structures is False and proteins_flag is False:
        weighted_bitscore_df = weighted_bitscore_df.with_columns(
            pl.col("query").str.split_exact(":", 1).alias(["contig_id", "cds_id"])
        ).drop(["query", "contig_id"])
    else:
        weighted_bitscore_df = weighted_bitscore_df.rename({"query": "cds_id"})

    # Drop columns based on conditions
    columns_to_drop = ["contig_id", "phrog", "product", "function"] if not proteins_flag else ["phrog", "product", "function"]
    filtered_tophits_df = filtered_tophits_df.drop(columns=columns_to_drop)

    # Merge dataframes
    merged_df = per_cds_df.join(filtered_tophits_df, on="cds_id", how="left")
    merged_df = merged_df.join(weighted_bitscore_df, on="cds_id", how="left")

    # Add annotation source and reorder columns
    product_index = merged_df.columns.index("product")
    column_order = [col for col in merged_df.columns[:product_index + 1] if col != "annotation_method"] + \
                ["annotation_method"] + \
                [col for col in merged_df.columns[product_index + 1:] if col != "annotation_method"]

    merged_df = merged_df.select(column_order)

    # Save the merged dataframe
    merged_df_path = Path(output) / f"{prefix}_per_cds_predictions.tsv"
    merged_df.write_csv(merged_df_path, sep="\t")

    # Create sub-database outputs
    sub_dbs_created = create_sub_db_outputs(merged_df, database, output)

    # If not proteins, generate and save function counts
    if not proteins_flag:
        contig_ids = merged_df.select("contig_id").unique().to_list()

        functions_list = []

        for contig in contig_ids:
            contig_df = merged_df.filter(pl.col("contig_id") == contig)
            cds_count = contig_df.height

            # Count functions
            function_counts = {
                "connector": contig_df.filter(pl.col("function") == "connector").height,
                "DNA, RNA and nucleotide metabolism": contig_df.filter(pl.col("function") == "DNA, RNA and nucleotide metabolism").height,
                "head and packaging": contig_df.filter(pl.col("function") == "head and packaging").height,
                "integration and excision": contig_df.filter(pl.col("function") == "integration and excision").height,
                "lysis": contig_df.filter(pl.col("function") == "lysis").height,
                "moron, auxiliary metabolic gene and host takeover": contig_df.filter(pl.col("function") == "moron, auxiliary metabolic gene and host takeover").height,
                "other": contig_df.filter(pl.col("function") == "other").height,
                "tail": contig_df.filter(pl.col("function") == "tail").height,
                "transcription regulation": contig_df.filter(pl.col("function") == "transcription regulation").height,
                "unknown function": contig_df.filter(pl.col("function") == "unknown function").height
            }

            # Count additional categories
            phrog_counts = {
                "VFDB_Virulence_Factors": contig_df.filter(pl.col("phrog") == "vfdb").height,
                "CARD_AMR": contig_df.filter(pl.col("phrog") == "card").height,
                "ACR_anti_crispr": contig_df.filter(pl.col("phrog") == "acr").height,
                "Defensefinder": contig_df.filter(pl.col("phrog") == "defensefinder").height,
                "Netflax": contig_df.filter(pl.col("phrog") == "netflax").height,
            }

            # Create the function count DataFrame
            count_list = [
                cds_count,
                *function_counts.values(),
                *phrog_counts.values()
            ]
            
            description_list = [
                "CDS",
                *function_counts.keys(),
                *phrog_counts.keys()
            ]

            contig_list = [contig] * len(count_list)

            cds_df = pl.DataFrame({
                "Description": description_list,
                "Count": count_list,
                "Contig": contig_list
            })

            functions_list.append(cds_df)

        description_total_df = pl.concat(functions_list)

        # Save the total description DataFrame
        descriptions_total_path = Path(output) / f"{prefix}_all_cds_functions.tsv"
        description_total_df.write_csv(descriptions_total_path, sep="\t")

    return sub_dbs_created