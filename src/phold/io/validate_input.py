from pathlib import Path

# imports
from loguru import logger

from phold.io.handle_genbank import get_fasta_run_pyrodigal_gv, get_genbank


def validate_input(input: Path, threads: int) -> dict:

    # validates input
    fasta_flag = False
    gb_dict = get_genbank(input)
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
        logger.info(f"Successfully parsed input {input} as a Genbank format file")

    return fasta_flag, gb_dict
