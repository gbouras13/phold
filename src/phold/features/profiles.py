import os
import shutil
from pathlib import Path

from Bio import SeqIO
from loguru import logger

from phold.features.create_foldseek_db import  foldseek_tsv2db
from phold.utils.util import remove_file


def generate_mmseqs_db_from_aa(
    fasta_aa: Path,  mmseqs2_db_path: Path, logdir: Path, prefix: str
) -> None:
    """
    Generate MMSeqs2 database from amino-acid sequences - for use with profiles later

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        mmseqs2_db_path (Path): Path to the directory where MMSeqs2 database will be stored.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.

    Returns:
        None
    """
    # read in amino-acid sequences
    sequences_aa = {}
    for record in SeqIO.parse(fasta_aa, "fasta"):
        sequences_aa[record.id] = str(record.seq)


    temp_aa_tsv: Path = Path(mmseqs2_db_path) / "aa.tsv"
    with open(temp_aa_tsv, "w") as f:
        for i,id in enumerate(sequences_aa.keys()):
            f.write("{}\t{}\n".format(str(i+1), sequences_aa[id]))

    temp_header_tsv: Path = Path(mmseqs2_db_path) / "header.tsv"
    with open(temp_header_tsv, "w") as f:
        for i,id in enumerate(sequences_aa.keys()):
            f.write("{}\t{}\n".format(str(i+1), id))

    # create MMSeqs2 db names
    short_db_name = f"{prefix}"
    aa_db_name: Path = Path(mmseqs2_db_path) / short_db_name
    header_db_name: Path = Path(mmseqs2_db_path) / f"{short_db_name}_h"

    # create Foldseek database with foldseek tsv2db

    foldseek_tsv2db(temp_aa_tsv, aa_db_name, 0, logdir)
    foldseek_tsv2db(temp_header_tsv, header_db_name, 12, logdir)

    # clean up
    remove_file(temp_aa_tsv)
    remove_file(temp_header_tsv)