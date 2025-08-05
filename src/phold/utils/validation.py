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
        logger.info(f"Successfully parsed input {input} as a {method} style Genbank file.")

    return fasta_flag, gb_dict, method


def instantiate_dirs(output_dir: Union[str, Path], force: bool) -> Path:
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
    try:
        process = sp.Popen(["foldseek", "version"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except:
        logger.error("Foldseek not found. Please reinstall phold.")

    foldseek_out, _ = process.communicate()
    foldseek_out = foldseek_out.decode()

    foldseek_version = foldseek_out.strip()

    if "941cd33" in foldseek_version:
        foldseek_major_version=10
        foldseek_minor_version="941cd33"
    else:
        logger.warning(f"Foldseek version found is v{foldseek_version}")
        logger.warning(f"Phold is recommended to be run with Foldseek v10.941cd33")
        logger.warning(f"Using a different Foldseek version is likely to work without issue, but this cannot be guaranteed.")

    logger.info(
        f"Foldseek version found is v{foldseek_major_version}.{foldseek_minor_version}"
    )

    #if foldseek_major_version != 10:
    #    logger.error("Foldseek is the wrong version. Please install v10.941cd33")
    #if foldseek_minor_version != "941cd33":
    #    logger.error("Foldseek is the wrong version. Please install v10.941cd33")

    logger.info("Foldseek version is ok")
