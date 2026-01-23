import os
import shutil
from pathlib import Path

from Bio import SeqIO
from loguru import logger

from phold.features.create_foldseek_db import  foldseek_tsv2db
from phold.utils.util import remove_file


def generate_mmseqs_db_from_aa(
    cds_dict: Path,  output: Path, logdir: Path, prefix: str, proteins_flag: bool
) -> None:
    """
    Generate MMSeqs2 database from amino-acid sequences - for use with profiles later

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        output (Path): Path to output directory.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.
        proteins_flag (bool): True if phold proteins-predict

    Returns:
        None
    """

    mmseqs2_db_path: Path = Path(output) / f"{prefix}_profiledb" # this will be the mmseqs2 db


    temp_aa_tsv = Path(output) / "aa.tsv"
    temp_header_tsv = Path(output) / "header.tsv"

    with open(temp_aa_tsv, "w") as aa_f, open(temp_header_tsv, "w") as h_f:
        idx = 1

        for contig_id, aa_contig_dict in cds_dict.items():
            if proteins_flag:
                for seq_id, aa_seq in aa_contig_dict.items():
                    aa_f.write(f"{idx}\t{aa_seq}\n")
                    h_f.write(f"{idx}\t{seq_id}\n")
                    idx += 1
            else:
                prefix = contig_id + ":"
                for seq_id, aa_seq in aa_contig_dict.items():
                    aa_f.write(f"{idx}\t{aa_seq}\n")
                    h_f.write(f"{idx}\t{prefix}{seq_id}\n")
                    idx += 1


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