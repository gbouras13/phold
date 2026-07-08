import re
import shutil
import subprocess as sp
import sys
from pathlib import Path
from typing import Dict, Union

from Bio import SeqIO
from Bio.SeqUtils import IUPACData
from loguru import logger

from phold.io.handle_genbank import get_fasta_run_pyrodigal_gv, get_genbank


def validate_input(input: Path, threads: int) -> Dict[str, Union[bool, Dict]]:
    """
    Validate the input file format and retrieve genomic data.

    Parameters:
        input (Path): Path to the input file.
        threads (int): Number of threads to use for prediction.

    Returns:
        Dict[str, Union[bool, Dict]]: A dictionary containing validation flags and genomic data.
    """

    # validates input
    fasta_flag = False
    gb_dict, method = get_genbank(input)
    if not gb_dict:
        logger.warning(f"{input} was not a Genbank format file")
        logger.warning(
            f"Now checking if the input {input} is a genome in nucleotide FASTA format"
        )
        logger.warning(f"pyrodigal-gv will be used to predict CDS")
        logger.warning(f"Phold will not predict tRNAs, tmRNAs or CRISPR repeats")
        logger.warning(
            f"Please use pharokka https://github.com/gbouras13/pharokka if you would like to predict these"
        )
        logger.warning(
            f"And then use the genbank output pharokka.gbk as --input for phold"
        )

        # check the contig ids are < 54 chars
        for record in SeqIO.parse(input, "fasta"):
            # Check if the length of the record ID is 54 characters or more
            if len(record.id) >= 54:
                logger.warning(
                    f"The contig header {record.id} is longer than 54 characters. It is recommended that you use shorter contig headers as this can create issues downstream."
                )

        gb_dict = get_fasta_run_pyrodigal_gv(input, threads)
        if not gb_dict:
            logger.warning("Error: no records found in FASTA file")
            logger.error("Please check your input")
        else:
            logger.info(
                f"Successfully parsed input {input} as a FASTA and predicted CDS"
            )
            fasta_flag = True
    else:
        logger.info(
            f"Successfully parsed input {input} as a {method} style Genbank file."
        )

    return fasta_flag, gb_dict, method


def instantiate_dirs(output_dir: Union[str, Path], force: bool, restart: bool = False) -> Path:
    """
    Checks and instantiates the output directory.

    Parameters:
        output_dir (Union[str, Path]): Path to the output directory.
        force (bool): Force flag indicating whether to overwrite existing directory.

    Returns:
        Path: Final output directory path.
    """

    # Checks the output directory
    # remove outdir on force
    logger.add(lambda _: sys.exit(1), level="ERROR")

    if restart and force:
        logger.warning(f"You have specified --restart and --force")
        logger.warning(f"This conflicts: proceeding assuming you want --restart and not --force")
        force = False


    if restart:
        logger.info(f"You have specified --restart")
        logger.info(f"Checking the output directory {output_dir} exists and contains foldseek_results.tsv")
        if Path(output_dir).exists() is False:
            logger.error(f"The output directory {output_dir} does not exist!")
        foldseek_results_path = Path(output_dir) / "foldseek_results.tsv"
        if Path(foldseek_results_path).exists() is False:
            logger.error(f"The foldseek_results.tsv file {foldseek_results_path} does not exist!")

        logger.info(f"The output directory {output_dir} exists and contains foldseek_results.tsv")

    else:

        logger.info(f"Checking the output directory {output_dir}")
        if force is True:
            if Path(output_dir).exists():
                logger.info(f"Removing {output_dir} because --force was specified")
                shutil.rmtree(output_dir)
            else:
                logger.info(
                    "--force was specified even though the output directory does not already exist. Continuing"
                )
        else:
            if Path(output_dir).exists():
                logger.error(
                    "Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory"
                )

        # instantiate outdir
        if Path(output_dir).exists() is False:
            Path(output_dir).mkdir(parents=True, exist_ok=True)


def check_dependencies() -> None:
    """
    Checks the dependencies and versions of non Python programs (i.e. Foldseek)

    Parameters:
        None

    Returns:
        None

    """

    #############
    # foldseek
    #############
    # Previously this used a bare ``except`` that logged "Foldseek not found"
    # and then fell through to ``process.communicate()`` on an undefined
    # ``process`` → the user's clear error message was followed by an
    # ``UnboundLocalError`` stack trace. Bare ``except`` also swallowed
    # Ctrl-C. Now: narrow the catch to the errors ``Popen`` actually raises
    # for a missing binary, surface the underlying error, and exit cleanly.
    try:
        process = sp.Popen(["foldseek", "version"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except (FileNotFoundError, PermissionError, OSError) as e:
        logger.error(
            f"Foldseek not found on PATH ({type(e).__name__}: {e}). "
            "Install foldseek (https://github.com/steineggerlab/foldseek) "
            "and ensure it is on your PATH, then re-run phold."
        )
        sys.exit(1)

    foldseek_out, _ = process.communicate()
    foldseek_out = foldseek_out.decode()

    foldseek_version = foldseek_out.strip()

    # Foldseek prints its version as ``<major>.<build_hash>`` (e.g.
    # ``10.941cd33``), sometimes just ``<major>``, occasionally with a ``v``
    # prefix. Phold is built against v10.x — any v10 build is expected to
    # work, but other majors aren't guaranteed. Parse the leading integer
    # as the major version; everything after is informational.
    EXPECTED_MAJOR = 10
    RECOMMENDED_VERSION = "10.941cd33"

    major_match = re.match(r"^v?(\d+)", foldseek_version)
    if major_match is None:
        # Unparseable: don't claim success, but don't block either —
        # phold will still try to run foldseek and the user will see a
        # more specific error if something is genuinely incompatible.
        logger.warning(
            f"Could not parse Foldseek version from {foldseek_version!r}. "
            f"Phold is built against Foldseek v{RECOMMENDED_VERSION}; "
            "continuing, but if you see Foldseek errors later this is "
            "the first thing to check."
        )
    elif int(major_match.group(1)) == EXPECTED_MAJOR:
        logger.info(
            f"Foldseek v{foldseek_version} detected (matches expected v{EXPECTED_MAJOR}.x)."
        )
    else:
        detected_major = int(major_match.group(1))
        logger.warning(
            f"Foldseek v{foldseek_version} detected (major version {detected_major}). "
            f"Phold is built and tested against Foldseek v{RECOMMENDED_VERSION}; "
            f"major version {detected_major} is likely to work but is not guaranteed."
        )
