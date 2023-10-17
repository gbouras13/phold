#!/usr/bin/env python3
"""

Adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/generate_foldseek_db.py

"""


from Bio import SeqIO
from pathlib import Path
from phold.utils.external_tools import ExternalTool
from phold.utils.util import remove_file
from loguru import logger


def generate_foldseek_db_from_aa_3di(fasta_aa: Path, fasta_3di: Path, foldseek_db_path: Path, logdir: Path, prefix: str ) -> None:


    # read in amino-acid sequences
    sequences_aa = {}
    for record in SeqIO.parse(fasta_aa, "fasta"):
        sequences_aa[record.id] = str(record.seq)

    # read in 3Di strings
    sequences_3di = {}
    for record in SeqIO.parse(fasta_3di, "fasta"):
        if not record.id in sequences_aa.keys():
            logger.warning("Warning: ignoring 3Di entry {}, since it is not in the amino-acid FASTA file".format(record.id))
        else:
            sequences_3di[record.id] = str(record.seq).upper()

    # assert that we parsed 3Di strings for all sequences in the amino-acid FASTA file
    for id in sequences_aa.keys():
        if not id in sequences_3di.keys():
            logger.warning("Warning: entry {} in amino-acid FASTA file has no corresponding 3Di string".format(id))
            logger.warning("Removing: entry {} from the Foldseek database ".format(id))
            sequences_aa = {id: sequence for id, sequence in sequences_aa.items() if id in sequences_3di}

    # generate TSV file contents
    tsv_aa = ""
    tsv_3di = ""
    tsv_header = ""
    for i,id in enumerate(sequences_aa.keys()):
        tsv_aa += "{}\t{}\n".format(str(i+1), sequences_aa[id])
        tsv_3di += "{}\t{}\n".format(str(i+1), sequences_3di[id])
        tsv_header += "{}\t{}\n".format(str(i+1), id)

    #### write temp tsv files

    # write TSV files
    temp_aa_tsv: Path = Path(foldseek_db_path) / "aa.tsv"
    temp_3di_tsv: Path = Path(foldseek_db_path) / "3di.tsv"
    temp_header_tsv: Path = Path(foldseek_db_path) / "header.tsv"
    with open(temp_aa_tsv, "w") as f:
        f.write(tsv_aa)
    with open(temp_3di_tsv, "w") as f:
        f.write(tsv_3di)
    with open(temp_header_tsv, "w") as f:
        f.write(tsv_header)

    
    # create foldseek db names

    short_db_name = f"{prefix}_foldseek_database"
    aa_db_name: Path = Path(foldseek_db_path) / short_db_name
    tsv_db_name: Path = Path(foldseek_db_path) / f"{short_db_name}_ss"
    header_db_name: Path = Path(foldseek_db_path) / f"{short_db_name}_h"

    # create Foldseek database with foldseek tsv2db

    foldseek_tsv2db(temp_aa_tsv, aa_db_name, 0, logdir) 
    foldseek_tsv2db(temp_3di_tsv, tsv_db_name, 0, logdir) 
    foldseek_tsv2db(temp_header_tsv, header_db_name, 12, logdir) 

    # clean up
    remove_file(temp_aa_tsv)
    remove_file(temp_3di_tsv)
    remove_file(temp_header_tsv)


def foldseek_tsv2db(in_tsv: Path, out_db_name: Path, db_type: int, logdir: Path) -> None:

    foldseek_tsv2db = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"tsv2db {in_tsv} {out_db_name}  --output-dbtype {str(db_type)} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_tsv2db)
