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

    mmseqs2_db_dir: Path = Path(output) / f"query_profiledb" # this will be subdir where the mmseqs2 db is 
    mmseqs2_db_dir.mkdir(parents=True, exist_ok=True)

    # create MMSeqs2 db names
    short_db_name = f"{prefix}"

    lookup_db_name: Path = Path(mmseqs2_db_dir) / f"{short_db_name}.lookup"

    temp_aa_tsv = Path(output) / "aa.tsv"
    temp_header_tsv = Path(output) / "header.tsv"

    with open(temp_aa_tsv, "w") as aa_f, open(temp_header_tsv, "w") as h_f, open(lookup_db_name, "w") as l_f:
        idx = 1

        for contig_id, aa_contig_dict in cds_dict.items():
            for seq_id, feature in aa_contig_dict.items():
                if proteins_flag:
                    header=seq_id
                else:
                    header = f"{contig_id}:{seq_id}"
                    aa_f.write(f"{idx}\t{feature.qualifiers["translation"]}\n")
                    h_f.write(f"{idx}\t{header}\n")
                    l_f.write(f"{idx-1}\t{header}\t0")
                    idx += 1

    
    aa_db_name: Path = Path(mmseqs2_db_dir) / short_db_name
    header_db_name: Path = Path(mmseqs2_db_dir) / f"{short_db_name}_h"

    # create Foldseek database with foldseek tsv2db

    foldseek_tsv2db(temp_aa_tsv, aa_db_name, 0, logdir)
    foldseek_tsv2db(temp_header_tsv, header_db_name, 12, logdir)

    # clean up
    remove_file(temp_aa_tsv)
    remove_file(temp_header_tsv)


def build_lookup(filename):
    """
    builds foldseek profile lookup using mmseqs lookup
    """
    lookup = {}
    with open(filename, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                key, value = parts[1], parts[0]
                lookup[key] = value
    return lookup
