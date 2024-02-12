#!/usr/bin/env python3


from datetime import datetime
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from loguru import logger

from phold.features.create_foldseek_db import generate_foldseek_db_from_aa_3di
from phold.features.query_remote_3Di import query_remote_3di
from phold.io.validate_input import validate_input
from phold.subcommands.compare import subcommand_compare
from phold.subcommands.predict import subcommand_predict
from phold.utils.util import (
    begin_phold,
    clean_up_temporary_files,
    end_phold,
    get_version,
    print_citation,
)
from phold.utils.validation import instantiate_dirs

# from phold.utils.validation import (
#     check_evalue,
#     instantiate_dirs,
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
            "-m",
            "--model_dir",
            required=False,
            type=click.Path(),
            default="model_dir",
            help="Path to save ProstT5_fp16 model to.",
        ),
        click.option(
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
        click.option(
            "--finetune",
            is_flag=True,
            help="Finetune",
        ),
        click.option(
            "--finetune_path", help="Path to finetuned model weights", default=None
        ),
        click.option(
            "--checkpoint_path", help="Path to CNN model weights", default=None
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
            help="Path to phold database.",
        ),
        click.option(
            "-e",
            "--evalue",
            default="1e-2",
            type=float,
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
            "--keep_tmp_files",
            is_flag=True,
            help="Keep temporary intermediate files, particularly the large foldseek_results.tsv of all Foldseek hits",
        ),
        click.option(
            "--split",
            is_flag=True,
            help="Split the Foldseek runs by ProstT5 probability",
        ),
        click.option(
            "--split_threshold",
            default=60,
            help="ProstT5 Probability to split by",
            type=float,
            show_default=True,
        ),
        click.option(
            "--card_vfdb_evalue",
            default="1e-10",
            type=float,
            help="stricter Evalue threshold for Foldseek CARD and VFDB hits",
            show_default=True,
        ),
        click.option(
            "--separate",
            is_flag=True,
            help="output separate genbank files for each contig",
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
    model_dir,
    model_name,
    database,
    batch_size,
    sensitivity,
    cpu,
    omit_probs,
    finetune,
    finetune_path,
    card_vfdb_evalue,
    split,
    split_threshold,
    checkpoint_path,
    separate,
    **kwargs,
):
    """phold predict then comapare all in one - GPU recommended"""

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
        "--batch_size": batch_size,
        "--sensitivity": sensitivity,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--finetune": finetune,
        "--finetune_path": finetune_path,
        "--checkpoint_path": checkpoint_path,
        "--split": split,
        "--split_threshold": split_threshold,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--separate": separate,
    }

    # initial logging etc
    start_time = begin_phold(params, "run")

    # validate input
    fasta_flag, gb_dict = validate_input(input, threads)

    # phold predict
    subcommand_predict(
        gb_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        batch_size,
        finetune,
        finetune_path,
        proteins_flag=False,
        checkpoint_path=checkpoint_path,
        fasta_flag=fasta_flag,
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
        pdb=False,
        pdb_dir=output,
        logdir=logdir,
        filter_pdbs=False,
        split=split,
        split_threshold=split_threshold,
        remote_flag=True,
        proteins_flag=False,
        separate=separate,
    )

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
    model_dir,
    model_name,
    batch_size,
    cpu,
    omit_probs,
    finetune,
    finetune_path,
    checkpoint_path,
    **kwargs,
):
    """Uses ProstT5 to predict 3Di tokens - GPU recommended"""

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
        "--batch_size": batch_size,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--finetune": finetune,
        "--finetune_path": finetune_path,
        "--checkpoint_path": checkpoint_path
    }

    # initial logging etc
    start_time = begin_phold(params, "predict")

    # validate input
    fasta_flag, gb_dict = validate_input(input, threads)

    # runs phold predict subcommand
    subcommand_predict(
        gb_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        batch_size,
        finetune,
        finetune_path,
        proteins_flag=False,
        checkpoint_path=checkpoint_path,
        fasta_flag=fasta_flag,
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
    "--pdb",
    is_flag=True,
    help="Use if you have pdbs for the input proteins (e.g. with AF2/Colabfold) that you specify with --pdb_dir",
)
@click.option(
    "--pdb_dir",
    help="Path to directory with pdbs.  The CDS IDs need to be in the name of the file",
    type=click.Path(),
)
@click.option(
    "--filter_pdbs",
    is_flag=True,
    help="Flag that creates a copy of the PDBs with matching record IDs found in the genbank. Helpful if you have a directory with lots of PDBs and want to annotate only e.g. 1 phage.",
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
    pdb,
    pdb_dir,
    filter_pdbs,
    keep_tmp_files,
    split,
    split_threshold,
    card_vfdb_evalue,
    separate,
    **kwargs,
):
    """Runs Foldseek vs phold db"""

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
        "--sensitivity": sensitivity,
        "--predictions_dir": predictions_dir,
        "--pdb": pdb,
        "--pdb_dir": pdb_dir,
        "--filter_pdbs": filter_pdbs,
        "--keep_tmp_files": keep_tmp_files,
        "--split": split,
        "--split_threshold": split_threshold,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--separate": separate,
    }

    # initial logging etc
    start_time = begin_phold(params, "compare")

    fasta_flag, gb_dict = validate_input(input, threads)

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
        pdb,
        pdb_dir,
        logdir,
        filter_pdbs,
        split,
        split_threshold,
        remote_flag=False,
        proteins_flag=False,
        fasta_flag=fasta_flag,
        separate=separate,
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
    model_dir,
    model_name,
    batch_size,
    cpu,
    omit_probs,
    finetune,
    finetune_path,
    checkpoint_path,
    **kwargs,
):
    """Runs ProstT5 on a multiFASTA input - GPU recommended"""

    # validates the directory  (need to before phold starts or else no log file is written)
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
        "--batch_size": batch_size,
        "--cpu": cpu,
        "--omit_probs": omit_probs,
        "--finetune": finetune,
        "--finetune_path": finetune_path,
        "--checkpoint_path": checkpoint_path
    }

    # initial logging etc
    start_time = begin_phold(params, "protein-predict")

    # validates fasta

    # Dictionary to store the records
    cds_dict = {}
    # need a dummmy nested dict
    cds_dict["proteins"] = {}

    # Iterate through the multifasta file and save each Seqfeature to the dictionary
    # 1 dummy record = proteins
    for record in SeqIO.parse(input, "fasta"):
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

    success = subcommand_predict(
        cds_dict,
        output,
        prefix,
        cpu,
        omit_probs,
        model_dir,
        model_name,
        batch_size,
        finetune,
        finetune_path,
        proteins_flag=True,
        checkpoint_path=checkpoint_path,
        fasta_flag=False,
    )

    # end phold
    end_phold(start_time, "proteins-predict")


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
    "--pdb",
    is_flag=True,
    help="Use if you have pdbs for the input proteins (e.g. with AF2/Colabfold) specified with --pdb_dir",
)
@click.option(
    "--pdb_dir",
    help="Path to directory with pdbs. The FASTA headers need to match names of the pdb files",
    type=click.Path(),
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
    pdb,
    pdb_dir,
    keep_tmp_files,
    split,
    split_threshold,
    card_vfdb_evalue,
    separate,
    **kwargs,
):
    """Runs Foldseek vs phold db on proteins input"""

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
        "--sensitivity": sensitivity,
        "--predictions_dir": predictions_dir,
        "--pdb": pdb,
        "--pdb_dir": pdb_dir,
        "--keep_tmp_files": keep_tmp_files,
        "--split": split,
        "--split_threshold": split_threshold,
        "--card_vfdb_evalue": card_vfdb_evalue,
    }

    # initial logging etc
    start_time = begin_phold(params, "proteins-compare")

    # validates fasta

    # Dictionary to store the records
    cds_dict = {}
    # need a dummmy nested dict
    cds_dict["proteins"] = {}

    # Iterate through the multifasta file and save each Seqfeature to the dictionary
    # 1 dummy record = proteins
    for record in SeqIO.parse(input, "fasta"):
        prot_id = record.id
        feature_location = FeatureLocation(0, len(record.seq))
        # Seq needs to be saved as the first element in list hence the closed brackets [str(record.seq)]
        seq_feature = SeqFeature(
            feature_location,
            type="CDS",
            qualifiers={
                "ID": record.id,
                "description": record.description,
                "translation": [str(record.seq)],
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
        pdb,
        pdb_dir,
        logdir,
        filter_pdbs=False,
        split=split,
        split_threshold=split_threshold,
        remote_flag=False,
        proteins_flag=True,
        fasta_flag=False,
        separate=False,
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
    split,
    split_threshold,
    card_vfdb_evalue,
    separate,
    **kwargs,
):
    """Uses foldseek API to run ProstT5 then foldseek locally"""

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
        "--sensitivity": sensitivity,
        "--keep_tmp_files": keep_tmp_files,
        "--split": split,
        "--split_threshold": split_threshold,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--separate": separate,
    }

    # initial logging etc
    start_time = begin_phold(params, "remote")

    # validate input
    fasta_flag, gb_dict = validate_input(input, threads)

    # Create a nested dictionary to store CDS features by contig ID
    cds_dict = {}

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"

    # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
    if fasta_flag is False:
        for record_id, record in gb_dict.items():
            cds_dict[record_id] = {}

            for cds_feature in record.features:
                if cds_feature.type == "CDS":
                    cds_dict[record_id][cds_feature.qualifiers["ID"][0]] = cds_feature
    else:  # if from pyrodigal, then exists
        cds_dict = gb_dict

    ## write the CDS to file

    with open(fasta_aa, "w+") as out_f:
        for contig_id, rest in cds_dict.items():
            aa_contig_dict = cds_dict[contig_id]

            # writes the CDS to file
            for seq_id, cds_feature in aa_contig_dict.items():
                out_f.write(f">{contig_id}:{seq_id}\n")
                out_f.write(f"{cds_feature.qualifiers['translation'][0]}\n")

    ############
    # prostt5 remote
    ############

    fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"
    query_remote_3di(cds_dict, fasta_3di)

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
        pdb=False,
        pdb_dir=output,
        logdir=logdir,
        filter_pdbs=False,
        split=split,
        split_threshold=split_threshold,
        remote_flag=True,
        proteins_flag=False,
        separate=separate,
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
    instantiate_dirs(output, force)

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
