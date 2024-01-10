#!/usr/bin/env python3
"""phold"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import click
from loguru import logger
import shutil
from Bio.SeqFeature import SeqFeature, FeatureLocation

from phold.io.handle_genbank import get_genbank, get_proteins

from phold.utils.util import begin_phold, end_phold, get_version, print_citation

from phold.features.predict_3Di import get_embeddings

from phold.features.create_foldseek_db import generate_foldseek_db_from_aa_3di, generate_foldseek_db_from_pdbs

from phold.results.tophit import get_tophits, parse_tophits, calculate_tophits_results

from phold.results.topfunction import (
    get_topfunctions,
    calculate_topfunctions_results,
)  # get_tophits, parse_tophits, calculate_tophits_results

from phold.results.distribution import get_distribution

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

"""
common options
"""

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
            help="Number of threads",
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
    ]
    for option in reversed(options):
        func = option(func)
    return func

"""
predict only options
"""

def predict_options(func):
    """predict command line args"""
    options = [
        click.option(
            "--model_dir",
            required=False,
            type=click.Path(),
            default="ProstT5_fp16_directory",
            help="Path to save ProstT5_fp16 model to.",
        ),
        click.option(
            "-m",
            "--model_name",
            required=False,
            type=str,
            default="Rostlab/ProstT5_fp16",
            show_default=True,
            help="Name of model: Rostlab/ProstT5_fp16.",
        ),
        click.option(
            "--batch_size",
            default=1,
            help="batch size for ProstT5. 1 is usually fastest.",
            show_default=True,
        ),
            click.option(
            "--cpu",
            is_flag=True,
            help="Use cpus only.",
        ),
            click.option(
            "--omit_probs",
            is_flag=True,
            help="Do not output 3Di probabilities from ProstT5",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func

"""
compare only options
"""

def compare_options(func):
    """compare command line args"""
    options = [
        click.option(
            "-d",
            "--database",
            required=True,
            type=click.Path(),
            help="Path to foldseek PHROGs or ENVHOGs database.",
        ),
        click.option(
            "--database_name",
            default="all_phrogs",
            type=str,
            required=False,
            show_default=True,
            help="Name of foldseek PHROGs or ENVHOGs database.",
        ),
        click.option(
            "-e",
            "--evalue",
            default="1e-3",
            help="e value threshold for Foldseek",
            show_default=True,
        ),
        click.option(
            "-s",
            "--sensitivity",
            default="9.5",
            help="sensitivity parameter for foldseek",
            type=float,
            show_default=True,
        ),
        click.option(
            "--mode",
            "mode",
            help="Mode to parse results.",
            default="tophit",
            show_default=True,
            type=click.Choice(["tophit", "topfunction", "distribution"]),
        ),
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
@predict_options
@compare_options
def run(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    model_dir,
    model_name,
    database,
    database_name,
    batch_size,
    sensitivity,
    mode,
    cpu,
    omit_probs,
    **kwargs,
):
    """Runs phold predict (ProstT5) and comapare (Foldseek)"""

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
        "--model_dir": model_dir,
        "--model_name": model_name,
        "--evalue": evalue,
        "--database": database,
        "--database_name": database_name,
        "--batch_size": batch_size,
        "--sensitivity": sensitivity,
        "--mode": mode,
        "--cpu": cpu,
        "--omit_probs": omit_probs
    }

    # initial logging etc
    start_time = begin_phold(params, "run")

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
            if cds_feature.type == "CDS":
                cds_dict[record_id][cds_feature.qualifiers["ID"][0]] = cds_feature

    ## write the CDS to file

    with open(fasta_aa, "w+") as out_f:
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

    if cpu is True:
        half_precision = False
    else:
        half_precision = True

    if omit_probs:
        output_probs = False
    else:
        output_probs = True

    get_embeddings(
        cds_dict,
        output,
        model_dir,
        model_name,
        half_precision=half_precision,
        max_residues=10000,
        max_seq_len=1000,
        max_batch=batch_size,
        proteins=False,
        cpu=cpu,
        output_probs=output_probs
    )

    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(
        fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
    )

    ###########
    # run foldseek search
    ###########

    short_db_name = prefix
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
        query_db, target_db, result_db, temp_db, threads, logdir, evalue, sensitivity
    )

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results.tsv"
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir)

    # tophits
    # calculate tophits

    if mode == "tophit":
        top_hits_tsv: Path = Path(output) / "tophits.tsv"
        filtered_tophits_df = get_tophits(result_tsv, top_hits_tsv, evalue)

        #####
        # parsing tophits depend on type of db
        #####

        filtered_tophits_df = parse_tophits(
            filtered_tophits_df, database, database_name
        )

        # calculate results and saves to tsvs

        calculate_tophits_results(filtered_tophits_df, cds_dict, output)

    elif mode == "topfunction":
        filtered_topfunctions_df = get_topfunctions(result_tsv, database, database_name)

        calculate_topfunctions_results(filtered_topfunctions_df, cds_dict, output)

    elif mode == "distribution":
        foldseek_df = get_distribution(result_tsv, database, database_name, output)

    # end phold
    end_phold(start_time, "run")


"""
predict command
Uses ProstT5 to predict 3Di sequences from AA, genbank
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
@predict_options
def predict(
    ctx,
    input,
    output,
    prefix,
    force,
    model_dir,
    model_name,
    batch_size,
    cpu,
    omit_probs,
    **kwargs,
):
    """Runs phold predict (ProstT5)"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--force": force,
        "--prefix": prefix,
        "--model_dir": model_dir,
        "--model_name": model_name,
        "--batch_size": batch_size,
        "--cpu": cpu,
        "--omit_probs": omit_probs
    }

    # initial logging etc
    start_time = begin_phold(params, "predict")

    # validates fasta
    gb_dict = get_genbank(input)
    if not gb_dict:
        logger.warning("Error: no sequences found in genbank file")
        logger.error("No sequences found in genbank file. Nothing to annotate")


    # Create a nested dictionary to store CDS features by contig ID
    cds_dict = {}
    fasta_aa: Path = Path(output) / "outputaa.fasta"

    # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
    for record_id, record in gb_dict.items():
        cds_dict[record_id] = {}

        for cds_feature in record.features:
            if cds_feature.type == "CDS":
                cds_dict[record_id][cds_feature.qualifiers["ID"][0]] = cds_feature
                

    ## write the CDS to file

    with open(fasta_aa, "w+") as out_f:
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
    # written to this file
    # output_3di: Path = Path(out_path) / "output3di.fasta"

    if cpu is True:
        half_precision = False
    else:
        half_precision = True

    if omit_probs:
        output_probs = False
    else:
        output_probs = True

    get_embeddings(
        cds_dict,
        output,
        model_dir,
        model_name,
        half_precision=half_precision,
        max_residues=10000,
        max_seq_len=1000,
        max_batch=batch_size,
        proteins=False,
        cpu=cpu,
        output_probs=output_probs
    )

    # end phold
    end_phold(start_time, "predict")


"""
compare command
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
@click.option(
    "--predictions_dir",
    help="Path to directory with output3di.faa",
    type=click.Path(),
)
@click.option(
    "--pdb",
    is_flag=True,
    help="Use if you have pdbs for the input proteins (with AF2/Colabfold).",
)
@click.option(
    "--pdb_dir",
    help="Path to directory with pdbs",
    type=click.Path(),
)
@click.option(
    "--unrelaxed",
    help="Use unrelaxed top rank pdb. By default, relaxed will be used (if found).",
    is_flag=True,
)
@common_options
@compare_options
def compare(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
    database_name,
    sensitivity,
    mode,
    predictions_dir,
    pdb,
    pdb_dir,
    unrelaxed,
    **kwargs,
):
    """Runs phold compare (Foldseek)"""

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
        "--database": database,
        "--database_name": database_name,
        "--sensitivity": sensitivity,
        "--mode": mode,
        "--predictions_dir": predictions_dir,
        "--pdb": pdb,
        "--pdb_dir": pdb_dir,
        "--unrelaxed": unrelaxed
    }

    # initial logging etc
    start_time = begin_phold(params, "compare")

    # validates fasta
    gb_dict = get_genbank(input)
    if not gb_dict:
        logger.warning("Error: no sequences found in genbank file")
        logger.error("No sequences found in genbank file. Nothing to annotate")

    # for key, value in gb_dict.items():
    #     logger.info(f"Parameter: {key} {value}.")

    # Create a nested dictionary to store CDS features by contig ID
    cds_dict = {}

    # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
    for record_id, record in gb_dict.items():
        cds_dict[record_id] = {}

        for cds_feature in record.features:
            if cds_feature.type == "CDS":
                cds_dict[record_id][cds_feature.qualifiers["ID"][0]] = cds_feature



    # # assumes these exist been run if pdb is false

    if pdb is False:
        fasta_aa_input: Path = Path(predictions_dir) / "outputaa.fasta"
        fasta_3di_input: Path = Path(predictions_dir) / "output3di.fasta"

    fasta_aa: Path = Path(output) / "outputaa.fasta"
    fasta_3di: Path = Path(output) / "output3di.fasta"


    ## copy the 3Di to file if pdb is false
    if pdb is False:
        if fasta_3di_input.exists():
            logger.info(f"Checked that the 3Di CDS file {fasta_3di_input} exists from phold predict.")
            shutil.copyfile(fasta_3di_input, fasta_3di)
        else:
            logger.error(f"The 3Di CDS file {fasta_3di} does not exist. Please run phold predict and/or check the prediction directory {output}.")
        # copy the aa to file 
        if fasta_aa_input.exists():
            logger.info(f"Checked that the AA CDS file {fasta_aa_input} exists from phold predict.")
            shutil.copyfile(fasta_aa_input, fasta_aa)
        else:
            logger.error(f"The AA CDS file {fasta_aa_input} does not exist. Please run phold predict and/or check the prediction directory {output}.")
    else:
        ## write the CDS to file
        logger.info(f"Writing the AAs to file {fasta_aa}.")
        with open(fasta_aa, "w+") as out_f:
            for contig_id, rest in cds_dict.items():
                aa_contig_dict = cds_dict[contig_id]

                # writes the CDS to file
                for seq_id, cds_feature in aa_contig_dict.items():
                    out_f.write(f">{contig_id}:{seq_id}\n")
                    out_f.write(f"{cds_feature.qualifiers['translation'][0]}\n")


    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    if pdb is True:
        logger.info("Creating a foldseek query db from the pdbs.")
        rank_001_pdb_path: Path = Path(output) / "rank_001_pdbs"
        rank_001_pdb_path.mkdir(parents=True, exist_ok=True)
        generate_foldseek_db_from_pdbs(fasta_aa, foldseek_query_db_path, pdb_dir, logdir, rank_001_pdb_path, prefix, unrelaxed
)
    else:
        generate_foldseek_db_from_aa_3di(
            fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
        )

    ###########
    # run foldseek search
    ###########

    short_db_name = prefix
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
        query_db, target_db, result_db, temp_db, threads, logdir, evalue, sensitivity
    )

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results.tsv"
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir)

    # tophits
    # calculate tophits

    if mode == "tophit":
        top_hits_tsv: Path = Path(output) / "tophits.tsv"
        filtered_tophits_df = get_tophits(result_tsv, top_hits_tsv, evalue)

        #####
        # parsing tophits depend on type of db
        #####

        filtered_tophits_df = parse_tophits(
            filtered_tophits_df, database, database_name
        )

        # calculate results and saves to tsvs
        calculate_tophits_results(filtered_tophits_df, cds_dict, output)

    elif mode == "topfunction":
        filtered_topfunctions_df = get_topfunctions(result_tsv, database, database_name)

        calculate_topfunctions_results(filtered_topfunctions_df, cds_dict, output)

    elif mode == "distribution":
        foldseek_df = get_distribution(result_tsv, database, database_name, output)

    # end phold
    end_phold(start_time, "run")



"""
proteins command
Uses ProstT5 to predict 3Di from a multiFASTA of proteins as input
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input multiFASTA file",
    type=click.Path(),
    required=True,
)
@common_options
@predict_options
def proteins(
    ctx,
    input,
    output,
    prefix,
    force,
    model_dir,
    model_name,
    batch_size,
    cpu,
    omit_probs,
    **kwargs,
):
    """Runs phold proteins (ProstT5 on a multiFASTA input)"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--force": force,
        "--prefix": prefix,
        "--model_dir": model_dir,
        "--model_name": model_name,
        "--batch_size": batch_size,
        "--cpu": cpu,
        "--omit_probs": omit_probs
    }

    # initial logging etc
    start_time = begin_phold(params, "protein")

    # validates fasta


    # Dictionary to store the records
    cds_dict = {}
    # need a dummmy nested dict
    cds_dict['proteins'] = {}


    # Iterate through the multifasta file and save each Seqfeature to the dictionary
    # 1 dummy record = proteins
    for record in SeqIO.parse(input, "fasta"):

        prot_id = record.id
        feature_location = FeatureLocation(0, len(record.seq))
        # Seq needs to be saved as the first element in list hence the closed brackets [str(record.seq)]
        seq_feature = SeqFeature(feature_location, 
                                 type="CDS", 
                                 strand=1,
                                 qualifiers={"ID": record.id,
                                             "description": record.description,
                                             "translation": [str(record.seq)]})
        
        
        cds_dict['proteins'][prot_id] = seq_feature

    if not cds_dict:
        logger.error(f"Error: no sequences found in {input} file")

    # copy input to output
    fasta_aa: Path = Path(output) / "outputaa.fasta"
    shutil.copyfile(input, fasta_aa)

    ############
    # prostt5
    ############

    # generates the embeddings using ProstT5 and saves them to file
    # written to this file
    output_3di: Path = Path(output) / "output3di.fasta"

    if cpu is True:
        half_precision = False
    else:
        half_precision = True

    if omit_probs:
        output_probs = False
    else:
        output_probs = True

    get_embeddings(
        cds_dict,
        output,
        model_dir,
        model_name,
        half_precision=half_precision,
        max_residues=10000,
        max_seq_len=1000,
        max_batch=batch_size,
        proteins=False,
        cpu=cpu,
        output_probs=output_probs
    )

    ## Need this to remove the fake record id 

    # Read the FASTA file
    records = list(SeqIO.parse(output_3di, "fasta"))

    # Process each record splitting the header to get rid of the fake record_id (DNE for this)
    for record in records:

        #header_parts = record.description.split(":")

        # Keep everything after the 10th character (to strip off proteins: )
        new_header = record.description[9:]
        #new_header = ":".join(header_parts[1:])

        # Update the record header
        record.id = new_header
        record.description = new_header

    # Write the modified records back to the same file
    with open(output_3di, "w") as output_file:
        SeqIO.write(records, output_file, "fasta")

    # end phold
    end_phold(start_time, "proteins")


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
@click.option(
    "-s",
    "--sensitivity",
    default="9.5",
    help="sensitivity parameter for foldseek",
    type=float,
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
    database,
    database_name,
    sensitivity,
    **kwargs,
):
    """Runs phold predict using foldseek API and compare locally"""

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
        "--database": database,
        "--database_name": database_name,
        "--sensitivity": sensitivity,
    }

    # initial logging etc
    start_time = begin_phold(params, "remote")

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
            if cds_feature.type == "CDS":
                cds_dict[record_id][cds_feature.qualifiers["ID"][0]] = cds_feature

    ## write the CDS to file

    with open(fasta_aa, "w+") as out_f:
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

    generate_foldseek_db_from_aa_3di(
        fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
    )

    ###########
    # run foldseek search
    ###########

    short_db_name = f"{prefix}_foldseek_database"
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
        query_db, target_db, result_db, temp_db, threads, logdir, evalue, sensitivity
    )

    # make result tsv
    result_tsv: Path = Path(output) / "foldseek_results.tsv"
    create_result_tsv(query_db, target_db, result_db, result_tsv, logdir)

    # calculate tophits
    top_hits_tsv: Path = Path(output) / "tophits.tsv"

    filtered_tophits_df = get_tophits(result_tsv, top_hits_tsv, evalue)

    # end phold
    end_phold(start_time, "remote")


"""
createphrog command 
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
    type=click.Path(),
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
    show_default=True,
)
@click.option(
    "--envhog_batch_size",
    type=int,
    default=1000,
    help="enhvog db protein batch size.",
    show_default=True,
)
def createphrog(
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
        "--envhog_batch_size": envhog_batch_size,
    }

    # initial logging etc
    start_time = begin_phold(params, "create")

    # gets proteins
    logger.info("Creating protein dictionary")
    prot_dict = get_proteins(input)
    logger.info("Protein dictionary created")

    if not prot_dict:
        logger.warning(f"Error: no sequences found in {input} FASTA file")
        logger.error(f"No sequences found in {input} FASTA file. Nothing to annotate")

    if envhog_flag is False:
        # check if the tsv exists
        tsv = Path(tsv)
        if tsv.exists() is False:
            logger.error(
                f"{tsv} does not exist. You need to specify a Path using --tsv to an input tsv linking Amino Acid FASTA file to phrog id"
            )

        # read tsv
        # needs to have 3 columns - seq_id, phrog and description
        phrog_mapping_df = pd.read_csv(
            tsv, sep="\t", names=["seq_id", "phrog", "description"]
        )

        # takes arguments gets all the desired phrogs
        phrog_list = ["phrog_" + str(i) for i in range(min_phrog, (max_phrog + 1))]

        # Create a nested dictionary to store CDS features by phrog ID
        cds_dict = {}

        logger.info("Creating dictionary to store all proteins for each phrog.")
        counter = 0

        # loops over all phrogs and adds them to dict - like a contig for normal phold
        for phrog_value in phrog_list:
            cds_dict[phrog_value] = {}

            # subset df
            phrog_group_df = phrog_mapping_df[phrog_mapping_df["phrog"] == phrog_value]
            # list of seq_ids
            seq_ids = phrog_group_df["seq_id"].tolist()

            # append to the cds dict
            for seq_id in seq_ids:
                if seq_id in prot_dict:
                    cds_dict[phrog_value][seq_id] = prot_dict[seq_id]
    else:
        envhog_end = envhog_start + envhog_batch_size - 1

        logger.info(
            f"You are running ProstT5 on enVhogs. Ignoring any PHROG specific commands."
        )
        logger.info(
            f"Taking the {envhog_start} to {envhog_end} proteins in your input file {input}."
        )
        # Convert the dictionary to a list of items (key-value pairs)
        dict_items = list(prot_dict.items())

        if envhog_end - 1 > len(dict_items):
            envhog_batch_size = len(dict_items) - envhog_end
            logger.warning(
                f"batch size reduced to {envhog_batch_size} to the end of the number of records."
            )
            envhog_end = len(dict_items)

        # Get items from 100th to 1000th
        entries = dict_items[envhog_start - 1 : envhog_end]

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
    with open(fasta_aa, "w+") as out_f:
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
    # required gpu
    fasta_3di: Path = Path(output) / "output3di.fasta"
    get_embeddings(
        cds_dict,
        output,
        model,
        half_precision=True,
        max_residues=3000,
        max_seq_len=500,
        max_batch=1,
        proteins=True,
        cpu=False
    )

    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output) / "foldseek_db"
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(
        fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
    )

    # end phold
    end_phold(start_time, "create")


"""
createdb command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "--fasta_aa",
    help="Path to input Amino Acid FASTA file of proteins",
    type=click.Path(),
    required=True,
)
@click.option(
    "--fasta_3di",
    help="Path to input 3Di FASTA file of proteins",
    type=click.Path(),
    required=True,
)
@click.option(
    "-o",
    "--output",
    default="output_phold_foldseek_db",
    show_default=True,
    type=click.Path(),
    help="Output directory ",
)
@click.option(
    "-t",
    "--threads",
    help="Number of threads to use with Foldseek",
    default=1,
    type=int,
    show_default=True,
)
@click.option(
    "-p",
    "--prefix",
    default="phold_foldseek_db",
    help="Prefix for Foldseek database",
    type=str,
    show_default=True,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Force overwrites the output directory",
)
def createdb(
    ctx,
    fasta_aa,
    fasta_3di,
    output,
    threads,
    prefix,
    force,
    **kwargs,
):
    """Creates phold compatible Foldseek db from AA FASTA and 3Di FASTA input files"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--aa_input": fasta_aa,
        "--fasta_3di": fasta_3di,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
    }

    # initial logging etc
    start_time = begin_phold(params, "createdb")

    logger.info(f"Creating the Foldseek database using {fasta_aa} and {fasta_3di}.")
    logger.info(
        f"The database will be saved in the {output} directory and be called {prefix}."
    )

    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output)
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    generate_foldseek_db_from_aa_3di(
        fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
    )

    # end phold
    end_phold(start_time, "createdb")


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
