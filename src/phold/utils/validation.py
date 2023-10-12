import re
import shutil
import subprocess as sp
import sys
from pathlib import Path

import click
from Bio import SeqIO
from loguru import logger


def instantiate_dirs(output_dir: Path, force: bool):
    """Checks the output directory
    :param out_dir: output directory path
    :param force: force flag
    :param logger: logger
    :return: out_dir: final output directory
    """

    # Checks the output directory
    # remove outdir on force
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Checking the output directory {output_dir}")
    if force is True:
        if Path(output_dir).exists():
            shutil.rmtree(output_dir)
        else:
            logger.info(
                "--force was specified even though the output directory does not already exist. Continuing."
            )
    else:
        if Path(output_dir).exists():
            logger.error(
                "Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory."
            )

    # instantiate outdir
    if Path(output_dir).exists() is False:
        Path(output_dir).mkdir(parents=True, exist_ok=True)



# def validate_fasta(input_fasta: Path):
#     """
#     Validates  FASTA input - that the input is a FASTA with 1 sequence
#     """
#     logger.info(
#         f"Checking that the input file {input_fasta} is in FASTA format and has only 1 entry."
#     )
#     # to get extension
#     with open(input_fasta, "r") as handle:
#         fasta = SeqIO.parse(handle, "fasta")
#         if any(fasta):
#             logger.info(f"{input_fasta} file checked.")
#         else:
#             logger.error(
#                 f"Error: {input_fasta} file is not in the FASTA format. Please check your input file"
#             )

#     with open(input_fasta, "r") as handle:
#         # Check the number of records
#         if len(list(SeqIO.parse(handle, "fasta"))) == 1:
#             logger.info(f"{input_fasta} has only one entry.")
#         else:
#             logger.error(
#                 f"{input_fasta} has more than one entry. Please check your input FASTA file!"
#             )



def is_protein_sequence(string):
    protein_letters = "acdefghiklmnpqrstvwy"
    nucleotide_letters = "acgnt"

    # Check if the string contains only nucleotide letters
    if all(letter.lower() in nucleotide_letters for letter in string):
        return False

    # Check if the string contains any protein letters
    return any(letter.lower() in protein_letters for letter in string)


def is_scientific_notation(evalue):
    """
    checks if evalue is scientific notation
    """
    # Define the regular expression pattern for scientific notation
    scientific_pattern = r"^[+\-]?(\d+(\.\d*)?|\.\d+)([eE][+\-]?\d+)?$"

    # Check if the number matches the scientific notation pattern
    return bool(re.match(scientific_pattern, evalue))


def is_numeric(evalue):
    """
    checks if evalue is numeric
    """
    try:
        float(evalue)  # Attempt to convert the value to a float
        return True
    except ValueError:
        return False


def check_evalue(evalue):
    """
    checks if the evalue is scientific notation or numeric
    """

    logger.info(f"You have specified an evalue of {evalue}.")

    if is_numeric(evalue) is False and is_scientific_notation(evalue) is False:
        logger.error(
            f"Error: evalue {evalue} is neither numeric nor in scientific notation. Please check your evalue."
        )

