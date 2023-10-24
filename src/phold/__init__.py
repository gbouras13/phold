#!/usr/bin/env python3
"""phold"""

import os
import shutil
from pathlib import Path
import pandas as pd
import copy

import click
from loguru import logger

from phold.io.handle_genbank import get_genbank, get_proteins

from phold.utils.util import (
    begin_phold,
    end_phold,
    get_version,
    print_citation
)

from phold.features.predict_3Di import get_embeddings

from phold.features.create_foldseek_db import  generate_foldseek_db_from_aa_3di

from phold.features.tophits import get_tophits

from phold.utils.validation import instantiate_dirs

from phold.features.run_foldseek import run_foldseek_search, create_result_tsv

from phold.features.query_remote_3Di import query_remote_3di

# from phold.utils.validation import (
#     check_evalue,
#     instantiate_dirs,
#     validate_choice_autocomplete,
#     validate_choice_mode,
#     validate_custom_db_fasta,
#     validate_fasta,
#     validate_fasta_all,
#     validate_fasta_bulk,
#     validate_ignore_file,
# )


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-o",
            "--output",
            default="output_phold",
            show_default=True,
            type=click.Path(),
            help="Output directory ",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads to use with Foldseek",
            default=1,
            type=int,
            show_default=True,
        ),
        click.option(
            "-p",
            "--prefix",
            default="phold",
            help="Prefix for output files",
            type=str,
            show_default=True,
        ),
        click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Force overwrites the output directory",
        ),
        click.option(
            "-m",
            "--model",
            required=False,
            type=str,
            default="Rostlab/ProstT5_fp16",
            help='Either a path to a directory holding the checkpoint for a pre-trained model or a huggingface repository link.' 
        ),
        click.option(
            "-d",
            "--database",
            required=True,
            type=click.Path(),
            help='Path to foldseek PHROGs database.' 
        )
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
def main_cli():
    1 + 1


"""
run command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
            "-i",
            "--input",
            help="Path to input file in Genbank format",
            type=click.Path(),
            required=True,
        )
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-3",
    help="e value threshold for Foldseek",
    show_default=True,
)
def run(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    model,
    database,
    **kwargs,
):
    """Runs phold"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--database": database
    }


    # initial logging etc
    start_time = begin_phold(params)

    # validates fasta
    gb_dict = get_genbank(input)
    if not gb_dict:
        logger.warning("Error: no sequences found in genbank file")
        logger.error("No sequences found in genbank file. Nothing to annotate")

    # for key, value in gb_dict.items():
    #     logger.info(f"Parameter: {key} {value}.")

    # Create a nested dictionary to store CDS features by contig ID
    cds_dict = {}

    fasta_aa: Path = Path(output) / "outputaa.fasta"

    # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
    
    for record_id, record in gb_dict.items():
        
        cds_dict[record_id] = {}

        for cds_feature in record.features:
            if cds_feature.type == 'CDS':
                cds_dict[record_id][cds_feature.qualifiers['ID'][0]] = cds_feature

    ## write the CDS to file

    with open(fasta_aa, 'w+') as out_f:
        for contig_id, rest in cds_dict.items():

            aa_contig_dict = cds_dict[contig_id]

            # writes the CDS to file
            for seq_id, cds_feature in aa_contig_dict.items():
                out_f.write(f">{contig_id}:{seq_id}\n")
                out_f.write(f"{cds_feature.qualifiers['translation'][0]}\n")

    ############
    # prostt5
    ############

    # generates the embeddings using ProstT5 and saves them to file
    fasta_3di: Path = Path(output) / "output3di.fasta"
    get_embeddings( cds_dict, output, model,  half_precision=True,    
                   max_residues=3000, max_seq_len=1000, max_batch=100 ) 
    
    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix )

    ###########
    # run foldseek search
    ###########

    short_db_name = f"{prefix}_foldseek_database"
    query_db: Path = Path(foldseek_query_db_path) / short_db_name
    target_db: Path = Path(database) / "all_phrog"

    # make result and temp dirs 
    result_db_base: Path = Path(output) / "result_db"
    result_db_base.mkdir(parents=True, exist_ok=True)
    result_db: Path = Path(result_db_base) / "result_db"

    temp_db: Path = Path(output) / "temp_db"
    temp_db.mkdir(parents=True, exist_ok=True)

    # run foldseek search
    run_foldseek_search(query_db, target_db,result_db, temp_db, threads, logdir, evalue )

    # make result tsv 
    result_tsv: Path =  Path(output) / "foldseek_results.tsv"
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir)

    # calculate tophits 
    top_hits_tsv: Path = Path(output) / "tophits.tsv"
    filtered_tophits_df = get_tophits(result_tsv, top_hits_tsv, evalue)
    # split the first column
    filtered_tophits_df[['record_id', 'cds_id']] = filtered_tophits_df['query'].str.split(':', expand=True, n=1)
    # Remove the original 'query' column
    filtered_tophits_df = filtered_tophits_df.drop(columns=['query'])

    # split on target
    filtered_tophits_df[['phrog', 'tophit_protein']] = filtered_tophits_df['target'].str.split(':', expand=True, n=1)
    # Remove the original 'target' column
    filtered_tophits_df = filtered_tophits_df.drop(columns=['target'])
    filtered_tophits_df['phrog'] = filtered_tophits_df['phrog'].str.replace('phrog_', '')
    filtered_tophits_df['phrog'] = filtered_tophits_df['phrog'].astype('str')

    # read in the mapping tsv

    phrog_annot_mapping_tsv: Path = Path(database) / "phrog_annot_v4.tsv"
    phrog_mapping_df = pd.read_csv(phrog_annot_mapping_tsv, sep='\t')
    phrog_mapping_df['phrog'] = phrog_mapping_df['phrog'].astype('str')

    # join the dfs

    filtered_tophits_df = filtered_tophits_df.merge(phrog_mapping_df, on='phrog', how='left')
    # Replace NaN values in the 'product' column with 'hypothetical protein'
    filtered_tophits_df['product'] = filtered_tophits_df['product'].fillna('hypothetical protein')

    # Convert the DataFrame to a nested dictionary
    result_dict = {}

    # instantiate the unique contig ids

    unique_contig_ids = filtered_tophits_df['record_id'].unique()

    for record_id in unique_contig_ids:
        result_dict[record_id] = {}
    
    for _, row in filtered_tophits_df.iterrows():
        record_id = row['record_id']
        cds_id = row['cds_id']
        values_dict = {
            'phrog': row['phrog'],
            'product': row['product'],
            'function': row['function'],
            'tophit_protein': row['tophit_protein'],
            'foldseek_alnScore': row['foldseek_alnScore'],
            'foldseek_seqIdentity': row['foldseek_seqIdentity'],
            'foldseek_eVal': row['foldseek_eVal'],
            'qStart': row['qStart'],
            'qEnd': row['qEnd'],
            'qLen': row['qLen'],
            'tStart': row['tStart'],
            'tEnd': row['tEnd'],
            'tLen': row['tLen'],
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

        original_functions_count_dict[record_id]['cds_count'] = len(updated_cds_dict[record_id])
        original_functions_count_dict[record_id]['phrog_count'] = 0
        original_functions_count_dict[record_id]['connector'] = 0
        original_functions_count_dict[record_id]['DNA, RNA and nucleotide metabolism'] = 0
        original_functions_count_dict[record_id]['head and packaging'] = 0
        original_functions_count_dict[record_id]['integration and excision'] = 0
        original_functions_count_dict[record_id]['lysis'] = 0
        original_functions_count_dict[record_id]['moron, auxiliary metabolic gene and host takeover'] = 0
        original_functions_count_dict[record_id]['other'] = 0
        original_functions_count_dict[record_id]['tail'] = 0
        original_functions_count_dict[record_id]['transcription regulation'] = 0
        original_functions_count_dict[record_id]['unknown function'] = 0

        new_functions_count_dict[record_id]['cds_count'] = len(updated_cds_dict[record_id])
        new_functions_count_dict[record_id]['phrog_count'] = 0
        new_functions_count_dict[record_id]['connector'] = 0
        new_functions_count_dict[record_id]['DNA, RNA and nucleotide metabolism'] = 0
        new_functions_count_dict[record_id]['head and packaging'] = 0
        new_functions_count_dict[record_id]['integration and excision'] = 0
        new_functions_count_dict[record_id]['lysis'] = 0
        new_functions_count_dict[record_id]['moron, auxiliary metabolic gene and host takeover'] = 0
        new_functions_count_dict[record_id]['other'] = 0
        new_functions_count_dict[record_id]['tail'] = 0
        new_functions_count_dict[record_id]['transcription regulation'] = 0
        new_functions_count_dict[record_id]['unknown function'] = 0
        new_functions_count_dict[record_id]['changed_phrogs'] = 0
        new_functions_count_dict[record_id]['same_phrogs'] = 0
        new_functions_count_dict[record_id]['foldseek_only_phrogs'] = 0
        new_functions_count_dict[record_id]['pharokka_only_phrogs'] = 0

        combined_functions_count_dict[record_id]['cds_count'] = len(updated_cds_dict[record_id])
        combined_functions_count_dict[record_id]['phrog_count'] = 0
        combined_functions_count_dict[record_id]['connector'] = 0
        combined_functions_count_dict[record_id]['DNA, RNA and nucleotide metabolism'] = 0
        combined_functions_count_dict[record_id]['head and packaging'] = 0
        combined_functions_count_dict[record_id]['integration and excision'] = 0
        combined_functions_count_dict[record_id]['lysis'] = 0
        combined_functions_count_dict[record_id]['moron, auxiliary metabolic gene and host takeover'] = 0
        combined_functions_count_dict[record_id]['other'] = 0
        combined_functions_count_dict[record_id]['tail'] = 0
        combined_functions_count_dict[record_id]['transcription regulation'] = 0
        combined_functions_count_dict[record_id]['unknown function'] = 0

        # iterates over the features
        # maybe can add 3DI as a genbank feature eventually?
        for cds_id, cds_feature in updated_cds_dict[record_id].items():
            # if pharokka got a phrog
            if cds_feature.qualifiers['phrog'][0] != "No_PHROG":
                original_functions_count_dict[record_id]['phrog_count'] += 1
            # get original function counts
            if cds_feature.qualifiers['function'][0] == "unknown function":
                original_functions_count_dict[record_id]['unknown function'] += 1
            elif cds_feature.qualifiers['function'][0] == "transcription regulation":
                original_functions_count_dict[record_id]['transcription regulation'] += 1
            elif cds_feature.qualifiers['function'][0] == "tail":
                original_functions_count_dict[record_id]['tail'] += 1
            elif cds_feature.qualifiers['function'][0] == "other":
                original_functions_count_dict[record_id]['other'] += 1      
            elif cds_feature.qualifiers['function'][0] == "moron":
                original_functions_count_dict[record_id]['moron, auxiliary metabolic gene and host takeover'] += 1 
            elif cds_feature.qualifiers['function'][0] == "lysis":
                original_functions_count_dict[record_id]['lysis'] += 1
            elif cds_feature.qualifiers['function'][0] == "integration and excision":
                original_functions_count_dict[record_id]['integration and excision'] += 1 
            elif cds_feature.qualifiers['function'][0] == "head and packaging":
                original_functions_count_dict[record_id]['head and packaging'] += 1
            elif cds_feature.qualifiers['function'][0] == "DNA":
                original_functions_count_dict[record_id]['DNA, RNA and nucleotide metabolism'] += 1
            elif cds_feature.qualifiers['function'][0] == "connector":
                original_functions_count_dict[record_id]['connector'] += 1 

            # now the updated dictionary
            if cds_id in result_dict[record_id].keys():
                # increase the phrog count
                new_functions_count_dict[record_id]['phrog_count'] += 1
                combined_functions_count_dict[record_id]['phrog_count'] += 1
                # update the counts
                if result_dict[record_id][cds_id]['function'] == "unknown function":
                    new_functions_count_dict[record_id]['unknown function'] += 1
                    combined_functions_count_dict[record_id]['unknown function'] += 1
                elif result_dict[record_id][cds_id]['function'] == "transcription regulation":
                    new_functions_count_dict[record_id]['transcription regulation'] += 1
                    combined_functions_count_dict[record_id]['transcription regulation'] += 1
                elif result_dict[record_id][cds_id]['function'] == "tail":
                    new_functions_count_dict[record_id]['tail'] += 1
                    combined_functions_count_dict[record_id]['tail'] += 1
                elif result_dict[record_id][cds_id]['function'] == "other":
                    new_functions_count_dict[record_id]['other'] += 1    
                    combined_functions_count_dict[record_id]['other'] += 1     
                elif result_dict[record_id][cds_id]['function'] == "moron, auxiliary metabolic gene and host takeover":
                    new_functions_count_dict[record_id]['moron, auxiliary metabolic gene and host takeover'] += 1 
                    combined_functions_count_dict[record_id]['moron, auxiliary metabolic gene and host takeover'] += 1 
                elif result_dict[record_id][cds_id]['function'] == "lysis":
                    new_functions_count_dict[record_id]['lysis'] += 1
                    combined_functions_count_dict[record_id]['lysis'] += 1
                elif result_dict[record_id][cds_id]['function'] == "integration and excision":
                    new_functions_count_dict[record_id]['integration and excision'] += 1 
                    combined_functions_count_dict[record_id]['integration and excision'] += 1 
                elif result_dict[record_id][cds_id]['function'] == "head and packaging":
                    new_functions_count_dict[record_id]['head and packaging'] += 1
                    combined_functions_count_dict[record_id]['head and packaging'] += 1
                elif result_dict[record_id][cds_id]['function'] == "DNA, RNA and nucleotide metabolism":
                    new_functions_count_dict[record_id]['DNA, RNA and nucleotide metabolism'] += 1
                    combined_functions_count_dict[record_id]['DNA, RNA and nucleotide metabolism'] += 1
                elif result_dict[record_id][cds_id]['function'] == "connector":
                    new_functions_count_dict[record_id]['connector'] += 1 
                    combined_functions_count_dict[record_id]['connector'] += 1 


                # update the phrog if different
                # same phrog
                if result_dict[record_id][cds_id]['phrog'] == cds_feature.qualifiers['phrog'][0]:
                    new_functions_count_dict[record_id]['same_phrogs'] += 1
                # different chrog
                if result_dict[record_id][cds_id]['phrog'] != cds_feature.qualifiers['phrog'][0]:
                    # where there was no phrog in pharokka
                    if cds_feature.qualifiers['phrog'][0] == "No_PHROG":
                        new_functions_count_dict[record_id]['foldseek_only_phrogs'] += 1
                    # different phrog to pharokka
                    else:
                        new_functions_count_dict[record_id]['changed_phrogs'] += 1
                    # update
                    updated_cds_dict[record_id][cds_id].qualifiers['phrog'][0] = result_dict[record_id][cds_id]['phrog']
                    updated_cds_dict[record_id][cds_id].qualifiers['product'][0] = result_dict[record_id][cds_id]['product']
                    updated_cds_dict[record_id][cds_id].qualifiers['function'][0] = result_dict[record_id][cds_id]['function']
            else: # will not be in results - unknown function
                new_functions_count_dict[record_id]['unknown function'] += 1
                if cds_feature.qualifiers['phrog'][0] != "No_PHROG":
                    new_functions_count_dict[record_id]['pharokka_only_phrogs'] +=1
                    combined_functions_count_dict[record_id]['phrog_count'] += 1
                    if cds_feature.qualifiers['function'][0] == "unknown function":
                        combined_functions_count_dict[record_id]['unknown function'] += 1
                    elif cds_feature.qualifiers['function'][0] == "transcription regulation":
                        combined_functions_count_dict[record_id]['transcription regulation'] += 1
                    elif cds_feature.qualifiers['function'][0] == "tail":
                        combined_functions_count_dict[record_id]['tail'] += 1
                    elif cds_feature.qualifiers['function'][0] == "other":
                        combined_functions_count_dict[record_id]['other'] += 1      
                    elif cds_feature.qualifiers['function'][0] == "moron":
                        combined_functions_count_dict[record_id]['moron, auxiliary metabolic gene and host takeover'] += 1 
                    elif cds_feature.qualifiers['function'][0] == "lysis":
                        combined_functions_count_dict[record_id]['lysis'] += 1
                    elif cds_feature.qualifiers['function'][0] == "integration and excision":
                        combined_functions_count_dict[record_id]['integration and excision'] += 1 
                    elif cds_feature.qualifiers['function'][0] == "head and packaging":
                        combined_functions_count_dict[record_id]['head and packaging'] += 1
                    elif cds_feature.qualifiers['function'][0] == "DNA":
                        combined_functions_count_dict[record_id]['DNA, RNA and nucleotide metabolism'] += 1
                    elif cds_feature.qualifiers['function'][0] == "connector":
                        combined_functions_count_dict[record_id]['connector'] += 1 
                else:
                    # no hits in either
                    combined_functions_count_dict[record_id]['unknown function'] += 1


                        
    print(original_functions_count_dict)
    print(new_functions_count_dict)
    print(combined_functions_count_dict)










    # validates fasta
    #validate_fasta(input)

    # validate e value
    #check_evalue(evalue)



    # end phold
    end_phold(start_time)


"""
remote command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
            "-i",
            "--input",
            help="Path to input file in Genbank format",
            type=click.Path(),
            required=True,
        )
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-3",
    help="e value threshold for Foldseek",
    show_default=True,
)
def remote(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    model,
    database,
    **kwargs,
):
    """Runs phold"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--database": database
    }


    # initial logging etc
    start_time = begin_phold(params)

    # validates fasta
    gb_dict = get_genbank(input)
    if not gb_dict:
        logger.warning("Error: no sequences found in genbank file")
        logger.error("No sequences found in genbank file. Nothing to annotate")

    # for key, value in gb_dict.items():
    #     logger.info(f"Parameter: {key} {value}.")

    # Create a nested dictionary to store CDS features by contig ID
    cds_dict = {}

    fasta_aa: Path = Path(output) / "outputaa.fasta"

    # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
    
    for record_id, record in gb_dict.items():
        
        cds_dict[record_id] = {}

        for cds_feature in record.features:
            if cds_feature.type == 'CDS':
                cds_dict[record_id][cds_feature.qualifiers['ID'][0]] = cds_feature

    ## write the CDS to file

    with open(fasta_aa, 'w+') as out_f:
        for contig_id, rest in cds_dict.items():

            aa_contig_dict = cds_dict[contig_id]

            # writes the CDS to file
            for seq_id, cds_feature in aa_contig_dict.items():
                out_f.write(f">{contig_id}:{seq_id}\n")
                out_f.write(f"{cds_feature.qualifiers['translation'][0]}\n")

    ############
    # prostt5
    ############

    # generates the embeddings using ProstT5 and saves them to file
    fasta_3di: Path = Path(output) / "output3di.fasta"
    
    query_remote_3di(cds_dict, output)
    
    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix )

    ###########
    # run foldseek search
    ###########

    short_db_name = f"{prefix}_foldseek_database"
    query_db: Path = Path(foldseek_query_db_path) / short_db_name
    target_db: Path = Path(database) #/ "toy_prophage_db"

    # make result and temp dirs 
    result_db_base: Path = Path(output) / "result_db"
    result_db_base.mkdir(parents=True, exist_ok=True)
    result_db: Path = Path(result_db_base) / "result_db"

    temp_db: Path = Path(output) / "temp_db"
    temp_db.mkdir(parents=True, exist_ok=True)

    # run foldseek search
    run_foldseek_search(query_db, target_db,result_db, temp_db, threads, logdir, evalue )

    # make result tsv 
    result_tsv: Path =  Path(output) / "foldseek_results.tsv"
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir)

    # calculate tophits 
    top_hits_tsv: Path = Path(output) / "tophits.tsv"

    filtered_tophits_df = get_tophits(result_tsv, top_hits_tsv, evalue)



    # validates fasta
    #validate_fasta(input)

    # validate e value
    #check_evalue(evalue)



    # end phold
    end_phold(start_time)


"""
create command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
            "-i",
            "--input",
            help="Path to input Amino Acid FASTA file in .faa",
            type=click.Path(),
            required=True,
        )
@click.option(
            "--tsv",
            help="Path to input tsv linking Amino Acid FASTA file to phrog id",
            default="phrog_mapping.tsv",
            type=click.Path()
        )
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-3",
    help="e value threshold for Foldseek",
    show_default=True,
)
@click.option(
    "--min_phrog",
    default=1,
    type=int,
    help="min phrog as integer (e.g. 1)",
    show_default=True,
)
@click.option(
    "--max_phrog",
    default=1000,
    type=int,
    help="max phrog as integer (e.g. 1000)",
    show_default=True,
)
@click.option(
    "--envhog_flag",
    is_flag=True,
    help="enhvog DB creation - assumes the input FASTA as already parsed.",
)
@click.option(
    "--envhog_start",
    type=int,
    default=1,
    help="enhvog db protein to start from.",
    show_default=True
)
@click.option(
    "--envhog_batch_size",
    type=int,
    default=1000,
    help="enhvog db protein batch size.",
    show_default=True
)

def create(
    ctx,
    input,
    output,
    tsv,
    threads,
    prefix,
    evalue,
    force,
    model,
    database,
    min_phrog,
    max_phrog,
    envhog_flag,
    envhog_start,
    envhog_batch_size,
    **kwargs,
):
    """Creates phold compatible PHROG or ENVHOG foldseek db using ProstT5"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--tsv": tsv,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--database": database,
        "--min_phrog": min_phrog,
        "--max_phrog": max_phrog,
        "--envhog_flag": envhog_flag,
        "--envhog_start": envhog_start,
        "--envhog_batch_size": envhog_batch_size

    }
    

    # initial logging etc
    start_time = begin_phold(params)

    # gets proteins
    logger.info('Creating protein dictionary')
    prot_dict = get_proteins(input)
    logger.info('Protein dictionary created')


    if not prot_dict:
        logger.warning(f"Error: no sequences found in {input} FASTA file")
        logger.error(f"No sequences found in {input} FASTA file. Nothing to annotate")

    if envhog_flag is False:

        # check if the tsv exists
        tsv = Path(tsv)
        if tsv.exists() is False:
            logger.error(f"{tsv} does not exist. You need to specify a Path using --tsv to an input tsv linking Amino Acid FASTA file to phrog id")

        # read tsv
        # needs to have 3 columns - seq_id, phrog and description
        phrog_mapping_df = pd.read_csv(tsv, sep='\t', names=["seq_id", "phrog", "description"])

        # takes arguments gets all the desired phrogs
        phrog_list =  ["phrog_" + str(i) for i in range(min_phrog, (max_phrog+1))]

        # Create a nested dictionary to store CDS features by phrog ID
        cds_dict = {}

        logger.info('Creating dictionary to store all proteins for each phrog.')
        counter = 0 


        # loops over all phrogs and adds them to dict - like a contig for normal phold
        for phrog_value in phrog_list:

            cds_dict[phrog_value] = {}

            # subset df
            phrog_group_df = phrog_mapping_df[phrog_mapping_df['phrog'] == phrog_value]
            # list of seq_ids
            seq_ids = phrog_group_df['seq_id'].tolist() 

            # append to the cds dict
            for seq_id in seq_ids:
                if seq_id in prot_dict:
                    cds_dict[phrog_value][seq_id] = prot_dict[seq_id]
    else:

        envhog_end = envhog_start + envhog_batch_size - 1

        logger.info(f"You are running ProstT5 on enVhogs. Ignoring any PHROG specific commands.")
        logger.info(f"Taking the {envhog_start} to {envhog_end} proteins in your input file {input}.")
        # Convert the dictionary to a list of items (key-value pairs)
        dict_items = list(prot_dict.items())

        if envhog_end - 1 > len(dict_items):
            envhog_batch_size = len(dict_items) - envhog_end
            logger.warning(f"batch size reduced to {envhog_batch_size} to the end of the number of records.")
            envhog_end = len(dict_items)

        # Get items from 100th to 1000th
        entries = dict_items[envhog_start-1:envhog_end]

        # to convert the selected items back to a dictionary:
        prot_dict = dict(entries)
        # need this hack to make it accept  the nested dictionary
        cds_dict = {}

        # get all unique phrogs in this set and instantiate the dict
        unique_phrogs = []
        for id, seq in prot_dict.items():
            parts = id.split(":")
            phrog = parts[0]
            if phrog not in unique_phrogs:
                unique_phrogs.append(phrog)

        for phrog in unique_phrogs:
            cds_dict[phrog] = {}

        for key, seq in prot_dict.items():
            parts = key.split(":")
            phrog = parts[0]
            prot_name = parts[1]
            cds_dict[phrog][prot_name] = prot_dict[key]

    
    fasta_aa: Path = Path(output) / "outputaa.fasta"

    ## write the CDS to file
    with open(fasta_aa, 'w+') as out_f:
        for phrog_id, rest in cds_dict.items():

            phrog_dict = cds_dict[phrog_id]

            # writes the CDS to file
            for seq_id, cds_feature in phrog_dict.items():
                out_f.write(f">{phrog_id}:{seq_id}\n")
                out_f.write(f"{phrog_dict[seq_id]}\n")




    ############
    # prostt5
    ############

    logger.info("Running ProstT5")

    # generates the embeddings using ProstT5 and saves them to file
    fasta_3di: Path = Path(output) / "output3di.fasta"
    get_embeddings( cds_dict, output, model,  half_precision=True,    
                   max_residues=3000, max_seq_len=500, max_batch=1, proteins=True ) 
    
    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix )

    # end phold
    end_phold(start_time)



@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


# main_cli.add_command(run)
main_cli.add_command(citation)


def main():
    main_cli()


if __name__ == "__main__":
    main()
