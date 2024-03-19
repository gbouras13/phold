#!/usr/bin/env python3
"""

Some code adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/generate_foldseek_db.py

"""


import os
import shutil
from pathlib import Path

from Bio import SeqIO
from loguru import logger

from phold.utils.external_tools import ExternalTool
from phold.utils.util import remove_file


def generate_foldseek_db_from_aa_3di(
    fasta_aa: Path, fasta_3di: Path, foldseek_db_path: Path, logdir: Path, prefix: str
) -> None:
    """
    Generate Foldseek database from amino-acid and 3Di sequences.

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        fasta_3di (Path): Path to the 3Di FASTA file.
        foldseek_db_path (Path): Path to the directory where Foldseek database will be stored.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.

    Returns:
        None
    """
    # read in amino-acid sequences
    sequences_aa = {}
    for record in SeqIO.parse(fasta_aa, "fasta"):
        sequences_aa[record.id] = str(record.seq)

    # read in 3Di strings
    sequences_3di = {}
    for record in SeqIO.parse(fasta_3di, "fasta"):
        if not record.id in sequences_aa.keys():
            logger.warning(
                "Warning: ignoring 3Di entry {}, since it is not in the amino-acid FASTA file".format(
                    record.id
                )
            )
        else:
            sequences_3di[record.id] = str(record.seq).upper()

    # assert that we parsed 3Di strings for all sequences in the amino-acid FASTA file
    for id in sequences_aa.keys():
        if not id in sequences_3di.keys():
            logger.warning(
                "Warning: entry {} in amino-acid FASTA file has no corresponding 3Di string".format(
                    id
                )
            )
            logger.warning("Removing: entry {} from the Foldseek database ".format(id))
            sequences_aa = {
                id: sequence
                for id, sequence in sequences_aa.items()
                if id in sequences_3di
            }

    # generate TSV file contents
    tsv_aa = ""
    tsv_3di = ""
    tsv_header = ""
    for i, id in enumerate(sequences_aa.keys()):
        tsv_aa += "{}\t{}\n".format(str(i + 1), sequences_aa[id])
        tsv_3di += "{}\t{}\n".format(str(i + 1), sequences_3di[id])
        tsv_header += "{}\t{}\n".format(str(i + 1), id)

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

    short_db_name = f"{prefix}"
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


def foldseek_tsv2db(
    in_tsv: Path, out_db_name: Path, db_type: int, logdir: Path
) -> None:
    """
    Convert a Foldseek TSV file to a Foldseek database.

    Args:
        in_tsv (Path): Path to the input TSV file.
        out_db_name (Path): Path for the output Foldseek database.
        db_type (int): Type of the output database.
        logdir (Path): Path to the directory where logs will be stored.

    Returns:
        None
    """
    foldseek_tsv2db = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"tsv2db {in_tsv} {out_db_name}  --output-dbtype {str(db_type)} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_tsv2db)


def generate_foldseek_db_from_pdbs(
    fasta_aa: Path,
    foldseek_db_path: Path,
    pdb_dir: Path,
    filtered_pdbs_path: Path,
    logdir: Path,
    prefix: str,
    filter_pdbs: bool,
) -> None:
    """
    Generate Foldseek database from PDB files.

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        foldseek_db_path (Path): Path to the directory where Foldseek database will be stored.
        pdb_dir (Path): Path to the directory containing PDB files.
        filtered_pdbs_path (Path): Path to the directory where filtered PDB files will be stored.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.
        filter_pdbs (bool): Flag indicating whether to filter PDB files or not.

    Returns:
        None
    """

    # read in amino-acid sequences
    sequences_aa = {}
    for record in SeqIO.parse(fasta_aa, "fasta"):
        sequences_aa[record.id] = str(record.seq)

    # lists all the pdb files
    pdb_files = [file for file in os.listdir(pdb_dir) if file.endswith(".pdb")]

    num_pdbs = len(pdb_files)

    num_pdbs = 0

    # Checks that ID is in the pdbs

    no_pdb_cds_ids = []

    for id in sequences_aa.keys():
        cds_id = id.split(":")[1]
        # record_id = id.split(":")[0]

        # this is potentially an issue if a contig has > 9999 AAs
        # need to fix with Pharokka possibly. Unlikely to occur but might!
        # enforce names as '{cds_id}.pdb'

        matching_files = [file for file in pdb_files if f"{cds_id}.pdb" == file]

        # delete the copying upon release, but for now do the copying to easy get the > Oct 2021 PDBs
        # with filter_pdbs
        if len(matching_files) == 1:
            if filter_pdbs is True:
                source_path = Path(pdb_dir) / matching_files[0]
                destination_path = Path(filtered_pdbs_path) / matching_files[0]
                shutil.copyfile(source_path, destination_path)
            num_pdbs += 1

        # should neve happen but in case
        if len(matching_files) > 1:
            logger.warning(f"More than 1 pdb found for {cds_id}")
            logger.warning("Taking the first one")
            if filter_pdbs is True:
                source_path = Path(pdb_dir) / matching_files[0]
                destination_path = Path(filtered_pdbs_path) / matching_files[0]
                shutil.copyfile(source_path, destination_path)
            num_pdbs += 1
        elif len(matching_files) == 0:
            logger.warning(f"No pdb found for {cds_id}")
            logger.warning(f"{cds_id} will be ignored in annotation")
            no_pdb_cds_ids.append(cds_id)

    if num_pdbs == 0:
        logger.error(
            f"No pdbs with matching CDS ids were found at all. Check the {pdb_dir} directory"
        )

    # generate the db
    short_db_name = f"{prefix}"
    pdb_db_name: Path = Path(foldseek_db_path) / short_db_name

    query_pdb_dir = pdb_dir

    # choose the filtered directory if true
    if filter_pdbs is True:
        query_pdb_dir = filtered_pdbs_path

    foldseek_createdb_from_pdbs = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createdb {query_pdb_dir} {pdb_db_name} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_createdb_from_pdbs)
