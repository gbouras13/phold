#!/usr/bin/env python3

from pathlib import Path

import click
from loguru import logger

# pycirclize is lazy-imported inside the `plot` subcommand (line ~1561) so
# the core phold pipeline doesn't drag pandas in transitively just to
# enable plotting. Users who run `phold plot` need `pip install pandas
# pycirclize` (pycirclize already requires pandas as a hard dep).

from phold.databases.db import (check_prostT5_download,
                                  download_zenodo_prostT5, install_database,
                                  validate_db)
from phold.features.create_foldseek_db import generate_foldseek_db_from_aa_3di
from phold.features.query_remote_3Di import query_remote_3di
# BioPython (SeqIO / SeqFeature) and phold.io.handle_genbank are lazy-imported
# in the two handler bodies that use them (proteins-predict and
# proteins-compare). handle_genbank transitively imports pyrodigal_gv.meta,
# which alone costs ~0.4 s on local disk and far more on a shared cluster
# filesystem — a price every subcommand, `--help` included, was paying.
# get_T5_model (from predict_3Di) and run_autotune (from autotune) are
# lazy-imported inside their handler bodies — both transitively import torch.
# subcommand_predict and subcommand_compare are lazy-imported inside each
# handler body (see below). subcommands.predict transitively imports torch
# (~4 s on a cold Python process), so importing it at module level made
# every phold subcommand — including install, citation, createdb, and plot —
# pay that cost even though they never call those functions.
from phold.utils.constants import CNN_DIR, DB_DIR
from phold.utils.util import (begin_phold, clean_up_temporary_files, end_phold,
                              get_version, print_citation)
from phold.utils.validation import (check_dependencies, instantiate_dirs,
                                    validate_input)
from importlib.resources import files

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

# ── model selection ───────────────────────────────────────────────────────────
# Kept as literal strings rather than importing pholdlib's registry at module
# level: every subcommand pays for imports here, and pholdlib.databases pulls in
# the torch-adjacent stack. resolve_model() imports lazily inside its body.
PROSTT5_MODEL = "prostt5"
MODERNPROST_MODEL_NAMES = [
    "modernprost-base",
    "modernprost-50M",
    "modernprost-pssm",
    "modernprost-50M-pssm",
]
MODEL_CHOICES = [PROSTT5_MODEL] + MODERNPROST_MODEL_NAMES


def resolve_model(model: str, database: Path, finetune: bool, vanilla: bool):
    """Turn the ``--model`` choice into the arguments subcommand_predict needs.

    Args:
        model: One of :data:`MODEL_CHOICES`.
        database: The validated Phold database directory, which doubles as the
            HuggingFace cache directory so model and DB live together.
        finetune: ProstT5 only — use the phage-finetuned encoder + CNN head.
        vanilla: ProstT5 only — use the CASP14-trained CNN head with the
            finetuned encoder.

    Returns:
        ``(model_name, model_dir, checkpoint_path, is_modernprost)``.
    """
    model = str(model).lower()
    # Choice(case_sensitive=False) lower-cases the value, so match back onto the
    # canonical registry spelling ("modernprost-50M", not "modernprost-50m").
    canonical = {name.lower(): name for name in MODERNPROST_MODEL_NAMES}

    if model in canonical:
        if finetune or vanilla:
            logger.warning(
                "--finetune and --vanilla only apply to --model prostt5; "
                f"ignoring them for {canonical[model]}."
            )
        return canonical[model], database, None, True

    model_name = "Rostlab/ProstT5_fp16"
    checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

    if finetune:
        model_name = "gbouras13/ProstT5Phold"
        checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt_finetune" / "phold_db_model.pth"
        if vanilla:
            checkpoint_path = (
                Path(CNN_DIR) / "cnn_chkpnt_finetune" / "vanilla_model.pth"
            )

    return model_name, database, checkpoint_path, False


def _resolved_task(model_name: str, task: str) -> str:
    """Resolve ``--task auto`` to the task *model_name* was trained for."""
    task = str(task).lower()
    if task != "auto":
        return task
    from pholdlib.databases.modernprost import default_task_for

    return default_task_for(model_name)

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
            "--model",
            "model",
            type=click.Choice(MODEL_CHOICES, case_sensitive=False),
            default="prostt5",
            show_default=True,
            help=(
                "Structure-token model. 'prostt5' is the ProstT5 encoder + CNN "
                "head (3Di only). The 'modernprost-*' models predict 3Di and a "
                "12-state alphabet in one pass and require a Phold database "
                "built with Foldseek 12-state support. The '-pssm' variants "
                "emit per-residue profiles and are searched as Foldseek "
                "profile databases."
            ),
        ),
        click.option(
            "--task",
            type=click.Choice(["auto", "classification", "pssm"], case_sensitive=False),
            default="auto",
            show_default=True,
            help=(
                "ModernProst inference mode. 'auto' uses whichever task the "
                "chosen model was trained for (pssm for the '-pssm' models, "
                "classification otherwise). Ignored for --model prostt5."
            ),
        ),
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
            "--max_batch_residues",
            type=int,
            default=None,
            help=(
                "Maximum residues per inference batch. Defaults to 5000 for "
                "ProstT5 and 50000 for ModernProst."
            ),
        ),
        click.option(
            "--cpu",
            is_flag=True,
            help="Use cpus only.",
        ),
        click.option(
            "--gpus",
            type=str,
            default=None,
            help=('Comma-separated CUDA device indices to use (e.g. "0,2"). '
                  "Default: all visible CUDA GPUs. Overridden by --cpu. "
                  "Has no effect on MPS / XPU systems."),
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
            default=0,
            help=(
                "Masks residues below this percentage model confidence before "
                "the Foldseek search. 0 disables masking. Applies to the 3Di "
                "and amino acid FASTAs for --model prostt5, and to the amino "
                "acid FASTA only for the ModernProst models, whose combined "
                "3Di+12-state alphabet has no masked state."
            ),
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
            "--evalue_12st_profile_comp",
            type=click.Choice(["0", "1", "2", "off"], case_sensitive=False),
            default="1",
            show_default=True,
            help=(
                "Composition source for Foldseek's 12-state e-value neural net "
                "on profile queries (--model modernprost-pssm / "
                "modernprost-50M-pssm). 1 reconstructs frequencies from the "
                "profile and uses those, rather than its centre sequence. "
                "0 = legacy whole-query, 2 = legacy whole-query 3Di+12-state, "
                "off = leave Foldseek's default. Only takes effect with "
                "Foldseek's NN e-value model (--evalue-nn-mode 2); no effect "
                "outside the profile path."
            ),
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
    model,
    task,
    autotune,
    batch_size,
    max_batch_residues,
    sensitivity,
    cpu,
    gpus,
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
    evalue_12st_profile_comp,
    custom_db,
    foldseek_gpu,
    hyps,
    finetune,
    vanilla,
    restart,
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
        "--model": model,
        "--task": task,
        "--autotune": autotune,
        "--batch_size": batch_size,
        "--max_batch_residues": max_batch_residues,
        "--sensitivity": sensitivity,
        "--keep_tmp_files": keep_tmp_files,
        "--cpu": cpu,
        "--gpus": gpus,
        "--omit_probs": omit_probs,
        "--card_vfdb_evalue": card_vfdb_evalue,
        "--separate": separate,
        "--max_seqs": max_seqs,
        "--save_per_residue_embeddings": save_per_residue_embeddings,
        "--save_per_protein_embeddings": save_per_protein_embeddings,
        "--ultra_sensitive": ultra_sensitive,
        "--mask_threshold": mask_threshold,
        "--extra_foldseek_params": extra_foldseek_params,
        "--evalue_12st_profile_comp": evalue_12st_profile_comp,
        "--custom_db": custom_db,
        "--foldseek_gpu": foldseek_gpu,
        "--hyps": hyps,
        "--finetune": finetune,
        "--vanilla": vanilla,
        "--restart": restart
    }

    from phold.subcommands.predict import subcommand_predict
    from phold.subcommands.compare import subcommand_compare
    from phold.features.autotune import run_autotune

    # initial logging etc
    start_time = begin_phold(params, "run")

    # check foldseek is installed
    check_dependencies()

    # check the database is installed and return it
    use_modernprost = model != PROSTT5_MODEL
    database = validate_db(database, DB_DIR, foldseek_gpu, require_12st=use_modernprost)

    # validate input
    fasta_flag, gb_dict, method = validate_input(input, threads)

    model_name, model_dir, checkpoint_path, use_modernprost = resolve_model(
        model, database, finetune, vanilla
    )
    profiles_flag = use_modernprost and _resolved_task(model_name, task) == "pssm"

    if not restart:
        # phold predict
        if autotune:
            if use_modernprost:
                # run_autotune builds a ProstT5 encoder + CNN head; it has no
                # ModernProst path. Skip rather than fail the whole run.
                logger.warning(
                    "--autotune is only supported for --model prostt5; "
                    f"using --batch_size {batch_size} for {model_name}."
                )
            else:
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
                    sample_seqs,
                    gpus=gpus,
                )

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
            gpus=gpus,
            task=task,
            max_batch_residues=max_batch_residues,
            logdir=logdir,
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
        evalue_12st_profile_comp=evalue_12st_profile_comp,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        restart=restart,
        gpus=gpus,
        ss_12st=use_modernprost,
        profiles=profiles_flag,
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
    model,
    task,
    autotune,
    batch_size,
    max_batch_residues,
    cpu,
    gpus,
    omit_probs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    mask_threshold,
    finetune,
    vanilla,
    hyps,
    **kwargs,
):
    """Predicts 3Di (and 12-state, for ModernProst models) tokens - GPU recommended"""

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
        "--model": model,
        "--task": task,
        "--autotune": autotune,
        "--batch_size": batch_size,
        "--max_batch_residues": max_batch_residues,
        "--cpu": cpu,
        "--gpus": gpus,
        "--omit_probs": omit_probs,
        "--save_per_residue_embeddings": save_per_residue_embeddings,
        "--save_per_protein_embeddings": save_per_protein_embeddings,
        "--mask_threshold": mask_threshold,
        "--finetune": finetune,
        "--vanilla": vanilla,
        "--hyps": hyps,
    }

    from phold.subcommands.predict import subcommand_predict
    from phold.features.autotune import run_autotune

    # initial logging etc
    start_time = begin_phold(params, "predict")

    # check the database is installed. The 12-state target DB is only needed by
    # `phold compare`, so predict alone does not require it — the model itself
    # is cached in the database directory either way.
    database = validate_db(database, DB_DIR, foldseek_gpu=False)

    # validate input
    fasta_flag, gb_dict, method = validate_input(input, threads)

    # runs phold predict subcommand
    model_name, model_dir, checkpoint_path, use_modernprost = resolve_model(
        model, database, finetune, vanilla
    )

    if autotune:
        if use_modernprost:
            logger.warning(
                "--autotune is only supported for --model prostt5; "
                f"using --batch_size {batch_size} for {model_name}."
            )
        else:
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
                sample_seqs,
                gpus=gpus,
            )

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
        gpus=gpus,
        task=task,
        max_batch_residues=max_batch_residues,
        logdir=logdir,
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
@click.option(
    "--gpus",
    type=str,
    default=None,
    help=('Comma-separated CUDA device indices for Foldseek-GPU (e.g. "0,2"). '
          "Default: all visible CUDA GPUs. Only meaningful with --foldseek_gpu."),
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
    gpus,
    keep_tmp_files,
    card_vfdb_evalue,
    separate,
    max_seqs,
    ultra_sensitive,
    extra_foldseek_params,
    evalue_12st_profile_comp,
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
        "--evalue_12st_profile_comp": evalue_12st_profile_comp,
        "--custom_db": custom_db,
        "--foldseek_gpu": foldseek_gpu,
        "--gpus": gpus,
        "--restart": restart
    }

    from phold.subcommands.compare import subcommand_compare

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
        evalue_12st_profile_comp=evalue_12st_profile_comp,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        restart=restart,
        gpus=gpus,
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
    model,
    task,
    autotune,
    batch_size,
    max_batch_residues,
    cpu,
    gpus,
    omit_probs,
    save_per_residue_embeddings,
    save_per_protein_embeddings,
    mask_threshold,
    finetune,
    vanilla,
    **kwargs,
):
    """Predicts 3Di (and 12-state, for ModernProst models) from a multiFASTA input - GPU recommended"""

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
        "--model": model,
        "--task": task,
        "--autotune": autotune,
        "--batch_size": batch_size,
        "--max_batch_residues": max_batch_residues,
        "--cpu": cpu,
        "--gpus": gpus,
        "--omit_probs": omit_probs,
        "--save_per_residue_embeddings": save_per_residue_embeddings,
        "--save_per_protein_embeddings": save_per_protein_embeddings,
        "--mask_threshold": mask_threshold,
        "--finetune": finetune,
        "--vanilla": vanilla,
    }

    from phold.subcommands.predict import subcommand_predict
    from phold.features.autotune import run_autotune

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

    from Bio import SeqIO
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    from phold.io.handle_genbank import open_protein_fasta_file

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
    model_name, model_dir, checkpoint_path, use_modernprost = resolve_model(
        model, database, finetune, vanilla
    )

    method = "pharokka"  # this can be whatever for proteins, it wont matter - it is for genbank input


    if autotune:
        if use_modernprost:
            logger.warning(
                "--autotune is only supported for --model prostt5; "
                f"using --batch_size {batch_size} for {model_name}."
            )
        else:
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
                sample_seqs,
                gpus=gpus,
            )

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
        gpus=gpus,
        task=task,
        max_batch_residues=max_batch_residues,
        logdir=logdir,
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
@click.option(
    "--gpus",
    type=str,
    default=None,
    help=('Comma-separated CUDA device indices for Foldseek-GPU (e.g. "0,2"). '
          "Default: all visible CUDA GPUs. Only meaningful with --foldseek_gpu."),
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
    gpus,
    keep_tmp_files,
    card_vfdb_evalue,
    max_seqs,
    ultra_sensitive,
    extra_foldseek_params,
    evalue_12st_profile_comp,
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
        "--evalue_12st_profile_comp": evalue_12st_profile_comp,
        "--custom_db": custom_db,
        "--foldseek_gpu": foldseek_gpu,
        "--gpus": gpus,
        "--restart": restart
    }

    from phold.subcommands.compare import subcommand_compare

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
    from Bio import SeqIO
    from Bio.SeqFeature import FeatureLocation, SeqFeature

    from phold.io.handle_genbank import open_protein_fasta_file

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
        evalue_12st_profile_comp=evalue_12st_profile_comp,
        custom_db=custom_db,
        foldseek_gpu=foldseek_gpu,
        restart=restart,
        gpus=gpus,
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
    evalue_12st_profile_comp,
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
        "--evalue_12st_profile_comp": evalue_12st_profile_comp,
        "--custom_db": custom_db,
    }

    from phold.subcommands.compare import subcommand_compare

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

    with open(fasta_aa, "w") as out_f:
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
        evalue_12st_profile_comp=evalue_12st_profile_comp,
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
    "--fasta_12st",
    help=(
        "Path to input 12-state FASTA file of proteins (e.g. phold predict's "
        "{prefix}_12st.fasta). When given, both alphabets are packed into one "
        "combined Foldseek database, which must be searched with --ss-12st 1. "
        "The 3Di FASTA must then be unmasked, as the combined encoding has no "
        "masked state."
    ),
    type=click.Path(),
    default=None,
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
    fasta_12st,
    output,
    threads,
    prefix,
    force,
    **kwargs,
):
    """Creates foldseek DB from AA FASTA and 3Di (and optionally 12-state) FASTA input files"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force, restart=False)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--fasta_aa": fasta_aa,
        "--fasta_3di": fasta_3di,
        "--fasta_12st": fasta_12st,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
    }

    # initial logging etc
    start_time = begin_phold(params, "createdb")

    # check foldseek is installed
    check_dependencies()

    if fasta_12st:
        logger.info(
            f"Creating the Foldseek database using {fasta_aa}, {fasta_3di} and "
            f"{fasta_12st}."
        )
        logger.info(
            "Both alphabets will be packed into one byte per residue "
            "(c = 3di_index * 12 + ss12_index). Search this database with "
            "--ss-12st 1."
        )
    else:
        logger.info(f"Creating the Foldseek database using {fasta_aa} and {fasta_3di}.")
    logger.info(
        f"The database will be saved in the {output} directory and be called {prefix}."
    )

    ############
    # create foldseek db
    ############

    foldseek_query_db_path: Path = Path(output)
    foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

    if fasta_12st:
        from phold.features.create_foldseek_db import \
            generate_foldseek_db_from_aa_3di_12st

        generate_foldseek_db_from_aa_3di_12st(
            fasta_aa, fasta_3di, fasta_12st, foldseek_query_db_path, logdir, prefix
        )
    else:
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
    "--model",
    "models",
    type=click.Choice(MODEL_CHOICES, case_sensitive=False),
    multiple=True,
    default=(PROSTT5_MODEL,),
    show_default=True,
    help=(
        "Model(s) to download. Repeat to install more than one, e.g. "
        "--model prostt5 --model modernprost-50M. The modernprost-base and "
        "modernprost-pssm checkpoints are ~4 GB each."
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
def install(
    ctx,
    database,
    foldseek_gpu,
    extended_db,
    models,
    threads,
    **kwargs,
):
    """Installs the structure-token model(s) and phold database"""

    from phold.features.predict_3Di import get_T5_model

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

    # always install with cpu mode as guarantee to be present
    cpu = True

    for model_choice in models:
        model_name, model_dir, _, use_modernprost = resolve_model(
            model_choice, database, finetune=False, vanilla=False
        )

        if use_modernprost:
            _install_modernprost(model_name, model_dir, threads)
            continue

        logger.info(
            f"Checking that the {model_name} ProstT5 model is available in {database}"
        )

    # Load (or download) the ProstT5 model. The check_fn / zenodo_fn
    # arguments are essential here: pholdlib's ``get_T5_model`` defaults
    # to ``local_files_only=True`` (offline) unless ``check_fn`` tells it
    # a download is needed. The first ``phold install`` on a fresh host
    # has nothing cached, so without check_fn the loader fails with
    # ``LocalEntryNotFoundError: outgoing traffic has been disabled``.
    # We pass the same pair that ``predict_3Di.get_embeddings`` uses, so
    # behaviour stays consistent across commands.
    #
    # ``get_T5_model`` returns ``(model, vocab, device)`` — the device
    # is irrelevant during install (we're only materialising the model
    # to disk), so we discard it with ``_``.
        model, vocab, _ = get_T5_model(
            database,
            model_name,
            cpu,
            threads=1,
            check_fn=check_prostT5_download,
            zenodo_fn=download_zenodo_prostT5,
        )
        del model
        del vocab
        logger.info(f"ProstT5 model downloaded")

    # will check if db is present, and if not, download it
    install_database(database, foldseek_gpu, extended_db, threads)


def _install_modernprost(model_name: str, model_dir: Path, threads: int) -> None:
    """Materialise a ModernProst checkpoint into the phold database directory.

    Split out of ``install`` so the torch-adjacent imports stay inside the
    function body — the same reason ``get_T5_model`` is imported lazily there.
    """
    from pholdlib.databases.modernprost import resolve_modernprost_model
    from pholdlib.modernprost.model import get_modernprost_model

    from phold.databases.db import (check_modernprost_download,
                                    modernprost_zenodo_downloader)

    spec = resolve_modernprost_model(model_name)
    logger.info(
        f"Checking that the {spec.hf_name} ModernProst model "
        f"(~{spec.download_mb} MB) is available in {model_dir}"
    )

    model, tokenizer, _ = get_modernprost_model(
        model_dir,
        model_name,
        cpu=True,
        threads=1,
        check_fn=check_modernprost_download,
        zenodo_fn=modernprost_zenodo_downloader(model_name),
    )
    del model
    del tokenizer
    logger.info(f"{spec.hf_name} model downloaded")


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

    # Lazy-imported here (not at module top) so the core phold pipeline
    # doesn't drag pandas in transitively just to support plotting.
    from pycirclize.parser import Genbank
    from phold.plot.plot import create_circos_plot

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
        if not Path(label_ids).exists():
            logger.error(f"{label_ids} was not found.")
        else:
            content = Path(label_ids).read_text()
            if not content.strip():
                logger.warning(f"{label_ids} contains no text. No CDS IDs will be force-labelled.")
            else:
                label_force_list = [
                    x.rstrip().split()[0] for x in content.splitlines() if x.strip()
                ]

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
    "--gpus",
    type=str,
    default=None,
    help=('Comma-separated CUDA device indices (e.g. "0,2"). '
          "Default: lowest visible CUDA GPU. Overridden by --cpu."),
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

def autotune(
    ctx,
    input,
    cpu,
    gpus,
    threads,
    database,
    step,
    min_batch,
    max_batch,
    sample_seqs,
    **kwargs,
):
    """Determines optimal batch size for 3Di prediction with your hardware"""

    params = {
        "--input": input,
        "--threads": threads,
        "--cpu": cpu,
        "--gpus": gpus,
        "--database": database,
        "--step": step,
        "--min_batch": min_batch,
        "--max_batch": max_batch,
        "--sample_seqs": sample_seqs,
    }

    from phold.features.autotune import run_autotune

    # initial logging etc
    start_time = begin_phold(params, "autotune")

    # check the database is installed
    database = validate_db(database, DB_DIR, foldseek_gpu=False)

    if input:
        input_path = input
    else:
        input_path = files("phold.features.autotune_data").joinpath("all_phold_structures_5000.fasta.gz")

    model_dir = database
    model_name = "Rostlab/ProstT5_fp16"

    batch_size = run_autotune(
        input_path,
        model_dir,
        model_name,
        cpu,
        threads,
        step,
        min_batch,
        max_batch,
        sample_seqs,
        gpus=gpus,
    )


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
