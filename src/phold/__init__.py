#!/usr/bin/env python3

import gzip
from pathlib import Path

import click
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from loguru import logger
from pycirclize.parser import Genbank

from phold.databases.db import install_database, validate_db
from phold.features.create_foldseek_db import generate_foldseek_db_from_aa_3di
from phold.features.predict_3Di import get_model
from phold.features.query_remote_3Di import query_remote_3di
from phold.io.handle_genbank import open_protein_fasta_file
from phold.plot.plot import create_circos_plot
from phold.subcommands.compare import subcommand_compare
from phold.subcommands.predict import subcommand_predict
from phold.utils.constants import CNN_DIR, DB_DIR
from phold.utils.util import (begin_phold, clean_up_temporary_files, end_phold,
                              get_version, print_citation)
from phold.utils.validation import (check_dependencies, instantiate_dirs,
                                    validate_input)
from phold.features.autotune import run_autotune
from importlib.resources import files

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
            "-d",
            "--database",
            type=str,
            default=None,
            help="Specific path to installed phold database",
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
            "--autotune",
            is_flag=True,
            help="Run autotuning to detect and automatically use best batch size for your hardware. Recommended only if you have a large dataset (e.g. thousands of proteins), or else autotuning will add rather than save runtime.",
        ),
        click.option(
            "--batch_size",
            default=1,
            help="batch size for ProstT5.",
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
            help="Do not output per residue 3Di probabilities from ProstT5. Mean per protein 3Di probabilities will always be output.",
        ),
        click.option(
            "--save_per_residue_embeddings",
            is_flag=True,
            help="Save the ProstT5 embeddings per resuide in a h5 file ",
        ),
        click.option(
            "--save_per_protein_embeddings",
            is_flag=True,
            help="Save the ProstT5 embeddings as means per protein in a h5 file",
        ),
        click.option(
            "--mask_threshold",
            default=25,
            help="Masks 3Di residues below this value of ProstT5 confidence for Foldseek searches",
            type=float,
            show_default=True,
        ),
        click.option(
            "--finetune",
            is_flag=True,
            help="Use gbouras13/ProstT5Phold encoder + CNN model both finetuned on phage proteins",
        ),
        click.option(
            "--vanilla",
            is_flag=True,
            help="Use vanilla CNN model (trained on CASP14) with ProstT5Phold encoder instead of the one trained on phage proteins",
        ),
        click.option(
            "--original",
            is_flag=True,
            help="Use original ProstT5 + CNN model, not ModernProst",
        ),
        click.option(
            "--hyps",
            is_flag=True,
            help="Use this to only annotate hypothetical proteins from a Pharokka GenBank input",
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
            "-e",
            "--evalue",
            default="1e-3",
            type=float,
            help="Evalue threshold for Foldseek",
            show_default=True,
        ),
        click.option(
            "-s",
            "--sensitivity",
            default="9.5",
            help="Sensitivity parameter for foldseek",
            type=float,
            show_default=True,
        ),
        click.option(
            "--keep_tmp_files",
            is_flag=True,
            help="Keep temporary intermediate files, particularly the large foldseek_results.tsv of all Foldseek hits",
        ),
        click.option(
            "--card_vfdb_evalue",
            default="1e-10",
            type=float,
            help="Stricter E-value threshold for Foldseek CARD and VFDB hits",
            show_default=True,
        ),
        click.option(
            "--separate",
            is_flag=True,
            help="Output separate GenBank files for each contig",
        ),
        click.option(
            "--max_seqs",
            type=int,
            default=1000,
            show_default=True,
            help="Maximum results per query sequence allowed to pass the prefilter. You may want to reduce this to save disk space for enormous datasets",
        ),
        click.option(
            "--ultra_sensitive",
            is_flag=True,
            help="Runs phold with maximum sensitivity by skipping Foldseek prefilter. Not recommended for large datasets.",
        ),
        click.option(
            "--extra_foldseek_params", type=str, help="Extra foldseek search params"
        ),
        click.option("--custom_db", type=str, help="Path to custom database"),
        click.option(
            "--foldseek_gpu",
            is_flag=True,
            help="Use this to enable compatibility with Foldseek-GPU search acceleration",
        ),
        click.option(
            "--restart",
            is_flag=True,
            help="Use this to restart phold from 'Processing Foldseek output' after foldseek_results.tsv is generated",
)
    ]
    for option in reversed(options):
        func = option(func)
    return func


"""
compare only options used for genbank/genome FASTA input (i.e. not proteins-compare)
"""


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
    help="Path to input file in Genbank format or nucleotide FASTA format",
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
    database,
    autotune,
    batch_size,
    sensitivity,
    cpu,
    omit_probs,
    keep_tmp_files,
    card_vfdb_evalue,
    separate,
    max_seqs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    ultra_sensitive,
    mask_threshold,
    extra_foldseek_params,
    custom_db,
    foldseek_gpu,
    hyps,
    finetune,
    vanilla,
    restart,
    original,
    **kwargs,
):
    """phold predict then comapare all in one - GPU recommended"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force, restart)

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
        "--autotune": autotune,
        "--batch_size": batch_size,
        "--sensitivity": sensitivity,
        "--keep_tmp_files": keep_tmp_files,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--separate": separate,
        "--max_seqs": max_seqs,
        "--save_per_residue_embeddings": save_per_residue_embeddings,
        "--save_per_protein_embeddings": save_per_protein_embeddings,
        "--ultra_sensitive": ultra_sensitive,
        "--mask_threshold": mask_threshold,
        "--extra_foldseek_params": extra_foldseek_params,
        "--custom_db": custom_db,
        "--foldseek_gpu": foldseek_gpu,
        "--hyps": hyps,
        "--finetune": finetune,
        "--vanilla": vanilla,
        "--restart": restart,
        "--original": original
    }

    # initial logging etc
    start_time = begin_phold(params, "run")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    database = validate_db(database, DB_DIR, foldseek_gpu)

    # validate input
    fasta_flag, gb_dict, method = validate_input(input, threads)


    if not restart:
        # phold predict

        model_dir = database
        model_name = "gbouras13/modernprost-base"
        checkpoint_path = None

        if original:
            model_name = "Rostlab/ProstT5_fp16"
            checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

        if finetune:
            model_name = "gbouras13/ProstT5Phold"
            checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "phold_db_model.pth"
            if vanilla:
                checkpoint_path = (
                    Path(CNN_DIR) / "cnn_chkpnt_finetune" / "vanilla_model.pth"
                )

        if autotune:

            input_path = files("phold.features.autotune_data").joinpath("all_phold_structures_5000.fasta.gz")

            step = 20
            min_batch = 1
            max_batch = 1001
            sample_seqs = 500

            batch_size = run_autotune(
                input_path,
                model_dir,
                model_name,
                cpu,
                threads,
                step, 
                min_batch,
                max_batch, 
                sample_seqs)

        subcommand_predict(
            gb_dict,
            method,
            output,
            prefix,
            cpu,
            omit_probs,
            model_dir,
            model_name,
            checkpoint_path,
            batch_size,
            proteins_flag=False,
            fasta_flag=fasta_flag,
            save_per_residue_embeddings=save_per_residue_embeddings,
            save_per_protein_embeddings=save_per_protein_embeddings,
            threads=threads,
            mask_threshold=mask_threshold,
            hyps=hyps,
        )

    # phold compare
    # predictions_dir is output as this will be where it lives
    subcommand_compare(
        gb_dict,
        output,
        threads,
        evalue,
        card_vfdb_evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir=output,
        structures=False,
        structure_dir=None,
        logdir=logdir,
        filter_structures=False,
        remote_flag=True,
        proteins_flag=False,
        fasta_flag=fasta_flag,
        separate=separate,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        restart=restart
    )

    # cleanup the temp files
    if keep_tmp_files is False:
        clean_up_temporary_files(output)

    # end phold
    end_phold(start_time, "run")


"""
predict command
Uses ProstT5 to predict 3Di sequences from AA, GenBank
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in Genbank format or nucleotide FASTA format",
    type=click.Path(),
    required=True,
)
@common_options
@predict_options
def predict(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    database,
    autotune,
    batch_size,
    cpu,
    omit_probs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    mask_threshold,
    finetune,
    vanilla,
    hyps,
    original,
    **kwargs,
):
    """Uses ProstT5 to predict 3Di tokens - GPU recommended"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force, restart=False)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--database": database,
        "--autotune": autotune,
        "--batch_size": batch_size,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--save_per_residue_embeddings": save_per_residue_embeddings,
        "--save_per_protein_embeddings": save_per_protein_embeddings,
        "--mask_threshold": mask_threshold,
        "--finetune": finetune,
        "--vanilla": vanilla,
        "--hyps": hyps,
        "--original": original
    }

    # initial logging etc
    start_time = begin_phold(params, "predict")

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu=False)

    # validate input
    fasta_flag, gb_dict, method = validate_input(input, threads)

    # runs phold predict subcommand
    model_dir = database

    model_name = "gbouras13/modernprost-base"
    checkpoint_path = None

    if original:
        model_name = "Rostlab/ProstT5_fp16"
        checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

    if finetune:
        model_name = "gbouras13/ProstT5Phold"
        checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "phold_db_model.pth"
        if vanilla:
            checkpoint_path = (
                Path(CNN_DIR) / "cnn_chkpnt_finetune" / "vanilla_model.pth"
            )

    if autotune:

        input_path = files("phold.features.autotune_data").joinpath("all_phold_structures_5000.fasta.gz")

        step = 20
        min_batch = 1
        max_batch = 1001
        sample_seqs = 500

        batch_size = run_autotune(
            input_path,
            model_dir,
            model_name,
            cpu,
            threads,
            step, 
            min_batch,
            max_batch, 
            sample_seqs)

    subcommand_predict(
        gb_dict,
        method,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        checkpoint_path,
        batch_size,
        proteins_flag=False,
        fasta_flag=fasta_flag,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
        hyps=hyps,
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
    help="Path to input file in Genbank format or nucleotide FASTA format",
    type=click.Path(),
    required=True,
)
@click.option(
    "--predictions_dir",
    help="Path to output directory from phold predict",
    type=click.Path(),
)
@click.option(
    "--structures",
    is_flag=True,
    help="Use if you have .pdb or .cif file structures for the input proteins (e.g. with AF2/Colabfold .pdb or AF3 for .cif) in a directory that you specify with --structure_dir",
)
@click.option(
    "--structure_dir",
    help="Path to directory with .pdb or .cif file structures. The CDS IDs need to be in the name of the file",
    type=click.Path(),
)
@click.option(
    "--filter_structures",
    is_flag=True,
    help="Flag that creates a copy of the .pdb or .cif files structures with matching record IDs found in the input GenBank file. Helpful if you have a directory with lots of .pdb files and want to annotate only e.g. 1 phage.",
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
    sensitivity,
    predictions_dir,
    structures,
    structure_dir,
    filter_structures,
    keep_tmp_files,
    card_vfdb_evalue,
    separate,
    max_seqs,
    ultra_sensitive,
    extra_foldseek_params,
    custom_db,
    foldseek_gpu,
    restart,
    **kwargs,
):
    """Runs Foldseek vs phold db"""

    # validates the directory  (need to before I start phold or else no log file is written)

    instantiate_dirs(output, force, restart)

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
        "--sensitivity": sensitivity,
        "--predictions_dir": predictions_dir,
        "--structures": structures,
        "--structure_dir": structure_dir,
        "--filter_structures": filter_structures,
        "--keep_tmp_files": keep_tmp_files,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--separate": separate,
        "--max_seqs": max_seqs,
        "--ultra_sensitive": ultra_sensitive,
        "--extra_foldseek_params": extra_foldseek_params,
        "--custom_db": custom_db,
        "--foldseek_gpu": foldseek_gpu,
        "--restart": restart
    }

    # initial logging etc
    start_time = begin_phold(params, "compare")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu)

    # validate fasta
    fasta_flag, gb_dict, method = validate_input(input, threads)

    subcommand_compare(
        gb_dict,
        output,
        threads,
        evalue,
        card_vfdb_evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir,
        structures,
        structure_dir,
        logdir,
        filter_structures,
        remote_flag=False,
        proteins_flag=False,
        fasta_flag=fasta_flag,
        separate=separate,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        restart=restart
    )

    # cleanup the temp files
    if keep_tmp_files is False:
        clean_up_temporary_files(output)

    # end phold
    end_phold(start_time, "compare")


""" 
proteins-predict command
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
def proteins_predict(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    database,
    autotune,
    batch_size,
    cpu,
    omit_probs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    mask_threshold,
    finetune,
    vanilla,
    original,
    **kwargs,
):
    """Runs ProstT5 on a multiFASTA input - GPU recommended"""

    # validates the directory  (need to before phold starts or else no log file is written)
    instantiate_dirs(output, force, restart=False)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--database": database,
        "--autotune": autotune,
        "--batch_size": batch_size,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--save_per_residue_embeddings": save_per_residue_embeddings,
        "--save_per_protein_embeddings": save_per_protein_embeddings,
        "--mask_threshold": mask_threshold,
        "--finetune": finetune,
        "--vanilla": vanilla,
        "--original": original
    }

    # initial logging etc
    start_time = begin_phold(params, "proteins-predict")

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu=False)

    # Dictionary to store the records
    cds_dict = {}
    # need a dummmy nested dict
    cds_dict["proteins"] = {}

    # Iterate through the multifasta file and save each Seqfeature to the dictionary
    # 1 dummy record = proteins

    with open_protein_fasta_file(input) as handle:  # handles gzip too
        records = list(SeqIO.parse(handle, "fasta"))
        if not records:
            logger.warning(f"No proteins were found in your input file {input}.")
            logger.error(
                f"Your input file {input} is likely not a amino acid FASTA file. Please check this."
            )
        for record in records:
            prot_id = record.id
            feature_location = FeatureLocation(0, len(record.seq))
            # Seq needs to be saved as the first element in list hence the closed brackets [str(record.seq)]
            seq_feature = SeqFeature(
                feature_location,
                type="CDS",
                qualifiers={
                    "ID": record.id,
                    "description": record.description,
                    "translation": str(record.seq),
                },
            )

            cds_dict["proteins"][prot_id] = seq_feature

    if not cds_dict:
        logger.error(f"Error: no AA protein sequences found in {input} file")

    # runs phold predict subcommand
    model_dir = database
    model_name = "gbouras13/modernprost-base"
    checkpoint_path = None

    if original:
        model_name = "Rostlab/ProstT5_fp16"
        checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"


    if finetune:
        model_name = "gbouras13/ProstT5Phold"
        checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "phold_db_model.pth"
        if vanilla:
            checkpoint_path = (
                Path(CNN_DIR) / "cnn_chkpnt_finetune" / "vanilla_model.pth"
            )

    method = "pharokka"  # this can be whatever for proteins, it wont matter - it is for genbank input


    if autotune:

        input_path = files("phold.features.autotune_data").joinpath("all_phold_structures_5000.fasta.gz")
        step = 20
        min_batch = 1
        max_batch = 1001
        sample_seqs = 500

        batch_size = run_autotune(
            input_path,
            model_dir,
            model_name,
            cpu,
            threads,
            step, 
            min_batch,
            max_batch, 
            sample_seqs)

    subcommand_predict(
        cds_dict,
        method,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        checkpoint_path,
        batch_size,
        proteins_flag=True,
        fasta_flag=False,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
        hyps=False,  # always False for this as no Pharokka genbank to parse on input
    )

    # end phold
    end_phold(start_time, "proteins-predict")


""" 
proteins compare command

Runs Foldseek vs phold DB for multiFASTA 3Di sequences (predicted with proteins-predict)
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in multiFASTA format",
    type=click.Path(),
    required=True,
)
@click.option(
    "--predictions_dir",
    help="Path to output directory from phold proteins-predict",
    type=click.Path(),
)
@click.option(
    "--structures",
    is_flag=True,
    help="Use if you have .pdb or .cif file structures for the input proteins (e.g. with AF2/Colabfold) in a directory that you specify with --structure_dir",
)
@click.option(
    "--structure_dir",
    help="Path to directory with .pdb or .cif file structures. The CDS IDs need to be in the name of the file",
    type=click.Path(),
)
@click.option(
    "--filter_structures",
    is_flag=True,
    help="Flag that creates a copy of the .pdb or .cif files structures with matching record IDs found in the input GenBank file. Helpful if you have a directory with lots of .pdb files and want to annotate only e.g. 1 phage.",
)
@common_options
@compare_options
def proteins_compare(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
    sensitivity,
    predictions_dir,
    structures,
    structure_dir,
    filter_structures,
    keep_tmp_files,
    card_vfdb_evalue,
    max_seqs,
    ultra_sensitive,
    extra_foldseek_params,
    custom_db,
    foldseek_gpu,
    restart,
    **kwargs
):
    """Runs Foldseek vs phold db on proteins input"""

    # validates the directory  (need to before I start phold or else no log file is written)

    instantiate_dirs(output, force, restart)

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
        "--sensitivity": sensitivity,
        "--predictions_dir": predictions_dir,
        "--structures": structures,
        "--structure_dir": structure_dir,
        "--filter_structures": filter_structures,
        "--keep_tmp_files": keep_tmp_files,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--max_seqs": max_seqs,
        "--ultra_sensitive": ultra_sensitive,
        "--extra_foldseek_params": extra_foldseek_params,
        "--custom_db": custom_db,
        "--foldseek_gpu": foldseek_gpu,
        "--restart": restart
    }

    # initial logging etc
    start_time = begin_phold(params, "proteins-compare")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu)

    # Dictionary to store the records
    cds_dict = {}
    # need a dummmy nested dict
    cds_dict["proteins"] = {}

    # Iterate through the multifasta file and save each Seqfeature to the dictionary
    # 1 dummy record = proteins
    with open_protein_fasta_file(input) as handle:  # handles gzip too
        records = list(SeqIO.parse(handle, "fasta"))
        if not records:
            logger.warning(f"No proteins were found in your input file {input}.")
            logger.error(
                f"Your input file {input} is likely not a amino acid FASTA file. Please check this."
            )
        for record in records:
            prot_id = record.id
            feature_location = FeatureLocation(0, len(record.seq))
            # Seq needs to be saved as the first element in list hence the closed brackets [str(record.seq)]
            seq_feature = SeqFeature(
                feature_location,
                type="CDS",
                qualifiers={
                    "ID": record.id,
                    "description": record.description,
                    "translation": str(record.seq),
                },
            )

            cds_dict["proteins"][prot_id] = seq_feature

    if not cds_dict:
        logger.error(f"Error: no AA protein sequences found in {input} file")

    success = subcommand_compare(
        cds_dict,
        output,
        threads,
        evalue,
        card_vfdb_evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir,
        structures,
        structure_dir,
        logdir,
        filter_structures,
        remote_flag=False,
        proteins_flag=True,
        fasta_flag=False,
        separate=False,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        restart=restart
    )

    # cleanup the temp files
    if keep_tmp_files is False:
        clean_up_temporary_files(output)

    # end phold
    end_phold(start_time, "proteins-compare")


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
    help="Path to input file in Genbank format or nucleotide FASTA format",
    type=click.Path(),
    required=True,
)
@common_options
@compare_options
def remote(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    database,
    sensitivity,
    keep_tmp_files,
    card_vfdb_evalue,
    separate,
    max_seqs,
    ultra_sensitive,
    extra_foldseek_params,
    custom_db,
    **kwargs,
):
    """Uses Foldseek API to run ProstT5 then Foldseek locally"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force, restart=False)

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
        "--sensitivity": sensitivity,
        "--keep_tmp_files": keep_tmp_files,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--separate": separate,
        "--max_seqs": max_seqs,
        "--ultra_sensitive": ultra_sensitive,
        "--extra_foldseek_params": extra_foldseek_params,
        "--custom_db": custom_db,
    }

    # initial logging etc
    start_time = begin_phold(params, "remote")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu=False)

    # validate input
    fasta_flag, gb_dict, method = validate_input(input, threads)

    # Create a nested dictionary to store CDS features by contig ID
    cds_dict = {}

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"

    # makes the nested dictionary {contig_id:{cds_id: cds_feature}}

    for record_id, record in gb_dict.items():
        cds_dict[record_id] = {}

        for cds_feature in record.features:
            if cds_feature.type == "CDS":
                if fasta_flag is False:
                    cds_feature.qualifiers["translation"] = cds_feature.qualifiers[
                        "translation"
                    ][0]
                    cds_dict[record_id][cds_feature.qualifiers["ID"][0]] = cds_feature
                else:
                    cds_dict[record_id][cds_feature.qualifiers["ID"]] = cds_feature

    ## write the CDS to file
    # FASTA -> takes the whole thing
    # Pharokka GBK -> requires just the first entry, the GBK is parsed as a list

    with open(fasta_aa, "w+") as out_f:
        for contig_id, rest in cds_dict.items():
            aa_contig_dict = cds_dict[contig_id]
            # writes the CDS to file
            for seq_id, cds_feature in aa_contig_dict.items():
                out_f.write(f">{contig_id}:{seq_id}\n")
                out_f.write(f"{cds_feature.qualifiers['translation']}\n")

    ############
    # prostt5 remote
    ############

    fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"
    query_remote_3di(cds_dict, fasta_3di, fasta_flag)

    ############
    # run compare vs db
    ############

    subcommand_compare(
        gb_dict,
        output,
        threads,
        evalue,
        card_vfdb_evalue,
        sensitivity,
        database,
        prefix,
        predictions_dir=output,
        structures=False,
        structure_dir=None,
        logdir=logdir,
        filter_structures=False,
        remote_flag=True,
        proteins_flag=False,
        fasta_flag=fasta_flag,
        separate=separate,
        max_seqs=max_seqs,
        ultra_sensitive=ultra_sensitive,
        extra_foldseek_params=extra_foldseek_params,
        custom_db=custom_db,
        foldseek_gpu=False,  # doesn't make sense for remote to do this as you wouldn't probably have a GPU
    )

    # cleanup the temp files
    if keep_tmp_files is False:
        clean_up_temporary_files(output)

    # end phold
    end_phold(start_time, "remote")


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
    """Creates foldseek DB from AA FASTA and 3Di FASTA input files"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force, restart=False)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--fasta_aa": fasta_aa,
        "--fasta_3di": fasta_3di,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
    }

    # initial logging etc
    start_time = begin_phold(params, "createdb")

    # check foldseek is installed
    check_dependencies()

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


"""
install command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-d",
    "--database",
    type=str,
    default=None,
    help="Specific path to install the phold database",
)
@click.option(
    "--foldseek_gpu",
    is_flag=True,
    help="Use this to enable compatibility with Foldseek-GPU acceleration",
)
@click.option(
    "--extended_db",
    is_flag=True,
    help=(
        "Download the extended Phold DB 3.16M including 1.8M efam and enVhog proteins without functional labels\n"
        "instead of the default Phold Search 1.36M. Using the extended database will likely marginally reduce\n"
        "functional annotation sensitivity and increase runtime, but may find more hits overall\n"
        "i.e. including to efam and enVhog proteins that have no functional labels."
    ),
)
@click.option(
    "-t",
    "--threads",
    help="Number of threads",
    default=1,
    type=int,
    show_default=True,
)
@click.option(
            "--original",
            is_flag=True,
            help="Use original ProstT5 + CNN model, not ModernProst",
        )
def install(
    ctx,
    database,
    foldseek_gpu,
    extended_db,
    threads,
    original,
    **kwargs,
):
    """Installs ProstT5 model and phold database"""

    if database is not None:
        logger.info(
            f"You have specified the {database} directory to store the Phold database and ProstT5 model"
        )
        database: Path = Path(database)
    else:
        logger.info(
            f"Downloading the Phold database into the default directory {DB_DIR}"
        )
        database = Path(DB_DIR)

    model_name = "gbouras13/modernprost-base"
    if original:
        model_name = "Rostlab/ProstT5_fp16"

    logger.info(
        f"Checking that the {model_name} model is available in {database}"
    )

    # always install with cpu mode as guarantee to be present
    cpu = True

    # load model (will be downloaded if not present)
    model, vocab = get_model(database, model_name, cpu, threads=1)
    del model
    del vocab
    logger.info(f"Model downloaded")

    # will check if db is present, and if not, download it
    install_database(database, foldseek_gpu, extended_db, threads)


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Path to input file in Genbank format (in the phold output directory)",
    type=click.Path(),
    required=True,
)
@click.option(
    "-o",
    "--output",
    default="phold_plots",
    show_default=True,
    type=click.Path(),
    help="Output directory to store phold plots",
)
@click.option(
    "-p",
    "--prefix",
    default="phold",
    help="Prefix for output files. Needs to match what phold was run with.",
    type=str,
    show_default=True,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Force overwrites the output directory",
)
@click.option("-a", "--all", is_flag=True, help="Plot every contig.")
@click.option(
    "-t",
    "--plot_title",
    default=None,
    help="Plot title. Only applies if --all is not specified. Will default to the phage's contig id.",
)
@click.option(
    "--label_hypotheticals",
    help="Flag to label hypothetical or unknown proteins. By default these are not labelled",
    is_flag=True,
)
@click.option(
    "--remove_other_features_labels",
    help="Flag to remove labels for tRNA/tmRNA/CRISPRs. By default these are labelled. \nThey will still be plotted in black",
    is_flag=True,
)
@click.option(
    "--title_size",
    type=float,
    default=20.0,
    help="Controls title size. Must be an integer. Defaults to 20",
)
@click.option(
    "--label_size",
    type=int,
    default=8,
    help="Controls annotation label size. Must be an integer. Defaults to 8",
)
@click.option(
    "--interval",
    default=5000,
    type=int,
    help="Axis tick interval. Must be an integer. Must be an integer. Defaults to 5000.",
)
@click.option(
    "--truncate",
    type=int,
    default=20,
    help="Number of characters to include in annoation labels before truncation with ellipsis. \nMust be an integer. Defaults to 20.",
)
@click.option(
    "--dpi",
    default="600",
    type=int,
    help="Resultion (dots per inch). Must be an integer. Defaults to 600.",
)
@click.option(
    "--annotations",
    default=1,
    type=float,
    help="Controls the proporition of annotations labelled. Must be a proportion between 0 and 1 inclusive. \n0 = no annotations, 0.5 = half of the annotations, 1 = all annotations. \nDefaults to 1. Chosen in order of CDS size.",
)
@click.option(
    "--label_ids",
    default=None,
    type=str,
    help="Text file with list of CDS IDs (from gff file) that are guaranteed to be labelled.",
)
def plot(
    ctx,
    prefix,
    input,
    output,
    force,
    all,
    plot_title,
    label_hypotheticals,
    remove_other_features_labels,
    title_size,
    label_size,
    interval,
    truncate,
    dpi,
    annotations,
    label_ids,
    **kwargs,
):
    """Creates Phold Circular Genome Plots"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--force": force,
        "--prefix": prefix,
        "--all": all,
        "--plot_title": plot_title,
        "--label_hypotheticals": label_hypotheticals,
        "--remove_other_features_labels": remove_other_features_labels,
        "--title_size": title_size,
        "--label_size": label_size,
        "--interval": interval,
        "--truncate": truncate,
        "--dpi": dpi,
        "--annotations": annotations,
        "--label_ids": label_ids,
    }

    # initial logging etc
    start_time = begin_phold(params, "plot")

    # single threaded plots
    threads = 1

    fasta_flag, gb_dict, method = validate_input(input, threads)

    gbk = Genbank(input)

    # gets all contigs and seqs
    gb_seq_dict = gbk.get_seqid2seq()

    gb_size_dict = gbk.get_seqid2size()

    contig_count = len(gb_seq_dict)

    # gets all features - will get all regardless of type (tRNA etc from pharokka)
    gb_feature_dict = gbk.get_seqid2features()

    # if there is 1 contig, then plot_title
    if contig_count > 1 and plot_title is not None:
        logger.warning(
            f"More than one contig found. Ignoring --plot_title {plot_title}"
        )

    # set contig id as title if single contig and no plot_title given
    if contig_count == 1 and plot_title is None:
        plot_title = str(contig_count)

    # check label_ids

    # list of all IDs that need to be labelled from file
    label_force_list = []

    if label_ids is not None:
        logger.info(
            f"You have specified a file {label_ids} containing a list of CDS IDs to force label."
        )
        # check if it is a file
        if Path(label_ids).exists() is False:
            logger.error(f"{label_ids} was not found.")
        # check if it contains text
        try:
            # Open the file in read mode
            with open(Path(label_ids), "r") as file:
                # Read the first character
                # will error if file is empty
                first_char = file.read(1)

                # read all the labels
                with open(Path(label_ids)) as f:
                    ignore_dict = {x.rstrip().split()[0] for x in f}
                # label force list
                label_force_list = list(ignore_dict)

        except FileNotFoundError:
            logger.warning(f"{label_ids} contains no text. No contigs will be ignored")

    # if there is 1 contig, then all the parameters will apply

    for contig_id, contig_sequence in gb_seq_dict.items():
        logger.info(f"Plotting {contig_id}")

        create_circos_plot(
            contig_id,
            contig_sequence,
            contig_count,
            gb_size_dict,
            gb_feature_dict,
            gbk,
            interval,
            annotations,
            title_size,
            plot_title,
            truncate,
            output,
            dpi,
            label_size,
            label_hypotheticals,
            remove_other_features_labels,
            label_force_list,
        )

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@click.option(
    "-i",
    "--input",
    help="Optional path to input file of proteins if you do not want to use the default sample of 5000 Phold DB proteins",
    type=click.Path()
)
@click.option(
    "--cpu",
    is_flag=True,
    help="Use cpus only.",
)
@click.option(
    "-t",
    "--threads",
    help="Number of threads",
    default=1,
    type=int,
    show_default=True,
)
@click.option(
    "-d",
    "--database",
    type=str,
    default=None,
    help="Specific path to installed phold database",
)
@click.option(
    "--min_batch",
    show_default=True,
    type=int,
    default=1,
    help="Minimum batch size to test",
)
@click.option(
    "--step",
    show_default=True,
    type=int,
    default=10,
    help="Controls batch size step increment",
)
@click.option(
    "--max_batch",
    default=251,
    show_default=True,
    type=int,
    help="Maximum batch size to test",
)
@click.option(
    "--sample_seqs",
    default=500,
    show_default=True,
    type=int,
    help="Number of proteins to subsample from input.",
)
@click.option(
            "--original",
            is_flag=True,
            help="Use original ProstT5 + CNN model, not ModernProst",
        )

def autotune(
    ctx,
    input,
    cpu,
    threads,
    database,
    step,
    min_batch,
    max_batch,
    sample_seqs,
    original,
    **kwargs,
):
    """Determines optimal batch size for 3Di prediction with your hardware"""

    params = {
        "--input": input,
        "--threads": threads,
        "--cpu": cpu,
        "--database": database,
        "--step": step,
        "--min_batch": min_batch,
        "--max_batch": max_batch,
        "--sample_seqs": sample_seqs,
        "--original": original
    }

    # initial logging etc
    start_time = begin_phold(params, "autotune")

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu=False)

    if input:
        input_path = input
    else:
        input_path = files("phold.features.autotune_data").joinpath("all_phold_structures_5000.fasta.gz")

    model_dir = database
    if original:
        model_name = "Rostlab/ProstT5_fp16"
    else:
         model_name = "gbouras13/modernprost-base"


    batch_size = run_autotune(
        input_path,
        model_dir,
        model_name,
        cpu,
        threads,
        step, 
        min_batch,
        max_batch, 
        sample_seqs)


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
