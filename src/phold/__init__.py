#!/usr/bin/env python3
"""phold"""

import os
import shutil
from pathlib import Path
import pandas as pd

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
    run_foldseek_search(query_db, target_db,result_db, temp_db, threads, logdir )

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
    target_db: Path = Path(database) / "toy_prophage_db"

    # make result and temp dirs 
    result_db_base: Path = Path(output) / "result_db"
    result_db_base.mkdir(parents=True, exist_ok=True)
    result_db: Path = Path(result_db_base) / "result_db"

    temp_db: Path = Path(output) / "temp_db"
    temp_db.mkdir(parents=True, exist_ok=True)

    # run foldseek search
    run_foldseek_search(query_db, target_db,result_db, temp_db, threads, logdir )

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
@click.option(
    "--min_phrog",
    required=True,
    type=int,
    help="min phrog as integer (e.g. 1)",
    show_default=True,
)
@click.option(
    "--max_phrog",
    required=True,
    type=int,
    help="max phrog as integer (e.g. 1000)",
    show_default=True,
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
    **kwargs,
):
    """Creates phold PHROG db using ProstT5"""

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
        "--max_phrog": max_phrog
    }


    # initial logging etc
    start_time = begin_phold(params)

    # gets proteins
    logger.info('Creating protein dictionary')
    prot_dict = get_proteins(input)
    logger.info('Protein dictionary created')


    if not prot_dict:
        logger.warning(f"Error: no sequences found in {input} FASTA file")
        logger.error("No sequences found in {input} FASTA file. Nothing to annotate")

    # read tsv
    # needs to have 3 columns - seq_id, phrog and description
    phrog_mapping_df = pd.read_csv(tsv, sep='\t', names=["seq_id", "phrog", "description"])

    # takes arguments gets all the desired phrogs
    phrog_list =  ["phrog_" + str(i) for i in range(min_phrog, (max_phrog+1))]

    # Create a nested dictionary to store CDS features by phrog ID
    cds_dict = {}

    logger.info('Creating dictionary to store all proteins for each phrog.')
    counter = 0 


    for phrog_value in phrog_list:

        counter += 1
        if counter % 100 == 0:
            print(counter)

        cds_dict[phrog_value] = {}

        # subset df
        phrog_group_df = phrog_mapping_df[phrog_mapping_df['phrog'] == phrog_value]
        # list of seq_ids
        seq_ids = phrog_group_df['seq_id'].tolist() 

        # append to the cds dict
        for seq_id in seq_ids:
            if seq_id in prot_dict:
                cds_dict[phrog_value][seq_id] = prot_dict[seq_id]


    
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
                   max_residues=3000, max_seq_len=1000, max_batch=100, proteins=True ) 
    
    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix )

    
    # validates fasta
    #validate_fasta(input)

    # validate e value
    #check_evalue(evalue)



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
