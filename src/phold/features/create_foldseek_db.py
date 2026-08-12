#!/usr/bin/env python3
"""

Some code adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/generate_foldseek_db.py

"""


import os
import shutil
from pathlib import Path
from typing import Dict

from Bio import SeqIO
from loguru import logger

# The pholdlib.modernprost builders are imported lazily inside the two
# functions that need them. phold/__init__.py imports this module at the top
# level for `generate_foldseek_db_from_aa_3di`, so an import here is paid by
# every subcommand — including `--help`. pholdlib's package __init__ is lazy
# too, but importing its submodules eagerly would still drag numpy in on that
# path, and the indirection is easy to lose track of; keeping the import next
# to its use makes the cost local and obvious.

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
            sequences_3di[record.id] = str(record.seq)  # no upper if masked

    # assert that we parsed 3Di strings for all sequences in the amino-acid FASTA file
    to_drop = set()
    for seq_id in sequences_aa.keys():
        if seq_id not in sequences_3di:
            logger.warning(
                "Warning: entry {} in amino-acid FASTA file has no corresponding 3Di string".format(
                    seq_id
                )
            )
            logger.warning("Removing: entry {} from the Foldseek database ".format(seq_id))
            to_drop.add(seq_id)
    for seq_id in to_drop:
        del sequences_aa[seq_id]

    # https://github.com/mheinzinger/ProstT5/issues/41

    temp_aa_tsv: Path = Path(foldseek_db_path) / "aa.tsv"
    with open(temp_aa_tsv, "w") as f:
        for i,id in enumerate(sequences_aa.keys()):
            f.write("{}\t{}\n".format(str(i+1), sequences_aa[id]))

    temp_3di_tsv: Path = Path(foldseek_db_path) / "3di.tsv"
    with open(temp_3di_tsv, "w") as f:
        for i,id in enumerate(sequences_aa.keys()):
            f.write("{}\t{}\n".format(str(i+1), sequences_3di[id]))

    temp_header_tsv: Path = Path(foldseek_db_path) / "header.tsv"
    with open(temp_header_tsv, "w") as f:
        for i,id in enumerate(sequences_aa.keys()):
            f.write("{}\t{}\n".format(str(i+1), id))
    #### write temp tsv files


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


def phold_tsv2db_runner(logdir: Path):
    """Build the ``tsv2db`` callback pholdlib's DB builders expect.

    pholdlib does not know about phold's ExternalTool wrapper, so it takes a
    ``(in_tsv, out_db, dbtype)`` callable. Binding *logdir* here keeps every
    ``foldseek tsv2db`` invocation logged in the same place as the rest of the
    run's external tool calls.
    """

    def _run(in_tsv: Path, out_db: Path, dbtype: int) -> None:
        foldseek_tsv2db(in_tsv, out_db, dbtype, logdir)

    return _run


def generate_foldseek_db_from_aa_3di_12st(
    fasta_aa: Path,
    fasta_3di: Path,
    fasta_12st: Path,
    foldseek_db_path: Path,
    logdir: Path,
    prefix: str,
) -> None:
    """
    Generate a combined 3Di + 12-state Foldseek database (ModernProst path).

    Both alphabets are packed into a single byte per residue in ``<prefix>_ss``
    (``c = 3di_index * 12 + ss12_index``). Search the result with
    ``--ss-12st 1``.

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        fasta_3di (Path): Path to the 3Di FASTA file (must be unmasked).
        fasta_12st (Path): Path to the 12-state FASTA file.
        foldseek_db_path (Path): Directory the Foldseek database is written to.
        logdir (Path): Directory where logs are stored.
        prefix (str): Prefix for the Foldseek database.

    Returns:
        None
    """
    from pholdlib.modernprost.foldseek_db import generate_combined_foldseek_db

    generate_combined_foldseek_db(
        fasta_aa,
        fasta_3di,
        fasta_12st,
        foldseek_db_path,
        prefix,
        phold_tsv2db_runner(logdir),
    )


def generate_foldseek_profile_db(
    profiles: Dict[str, Dict[str, Dict]],
    fasta_aa: Path,
    profile_db_path: Path,
    logdir: Path,
    prefix: str,
) -> Path:
    """
    Generate Foldseek 3Di / 12-state / amino-acid profile databases.

    Used by the ``-pssm`` ModernProst checkpoints, whose per-residue softmax
    distributions are scored directly rather than collapsed to a single state.

    Args:
        profiles (Dict): Nested ``{contig_id: {seq_id: {"3di": ndarray,
            "12st": ndarray}}}`` from ``get_modernprost_predictions``.
        fasta_aa (Path): Path to the amino-acid FASTA written by phold predict.
            Supplies both the sequences and the header formatting, so the
            profile DB keys match what the combined-DB path would produce.
        profile_db_path (Path): Directory the profile databases are written to.
        logdir (Path): Directory where logs are stored.
        prefix (str): Prefix for the Foldseek database.

    Returns:
        Path: prefix of the query profile DB to hand to ``foldseek search``.
    """
    from pholdlib.modernprost.foldseek_db import (generate_sequence_foldseek_db,
                                                  read_fasta)
    from pholdlib.modernprost.profile_db import write_profile_foldseek_dbs

    aa_sequences = read_fasta(fasta_aa)

    # Re-key the nested profiles onto the FASTA headers ("contig:cds" or plain
    # "cds"), which is what the amino-acid FASTA — and therefore the DB lookup
    # — uses. Anything the FASTA does not contain was dropped upstream (failed
    # inference, zero-length prediction) and must be dropped here too.
    flat_3di: Dict = {}
    flat_12st: Dict = {}
    for contig_id, contig_profiles in profiles.items():
        for seq_id, heads in contig_profiles.items():
            header = seq_id if seq_id in aa_sequences else f"{contig_id}:{seq_id}"
            if header not in aa_sequences:
                logger.warning(
                    f"Skipping {seq_id} in the Foldseek profile database: it has "
                    f"no entry in {fasta_aa}"
                )
                continue
            flat_3di[header] = heads["3di"]
            flat_12st[header] = heads["12st"]

    # Only build the DB over records that have a profile, so the amino-acid
    # profile DB cannot end up with keys the structural profile DBs lack.
    aa_sequences = {k: v for k, v in aa_sequences.items() if k in flat_3di}
    if not aa_sequences:
        logger.error(
            "No ModernProst profiles could be matched to the amino-acid FASTA; "
            "cannot build a Foldseek profile database."
        )

    profile_db_path = Path(profile_db_path)
    profile_db_path.mkdir(parents=True, exist_ok=True)

    source_db, lookup = generate_sequence_foldseek_db(
        aa_sequences, profile_db_path, prefix, phold_tsv2db_runner(logdir)
    )

    write_profile_foldseek_dbs(
        flat_3di,
        flat_12st,
        aa_sequences,
        lookup,
        profile_db_path,
        prefix,
        source_db,
    )

    return profile_db_path / f"{prefix}_profile"


def generate_foldseek_db_from_structures(
    fasta_aa: Path,
    foldseek_db_path: Path,
    structure_dir: Path,
    filtered_structures_path: Path,
    logdir: Path,
    prefix: str,
    filter_structures: bool,
    proteins_flag: bool,
) -> None:
    """
    Generate Foldseek database from PDB files.

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        foldseek_db_path (Path): Path to the directory where Foldseek database will be stored.
        structure_dir (Path): Path to the directory containing .pdb or .cif structure files.
        filtered_structures_path (Path): Path to the directory where filtered .pdb or .cif structure files will be stored.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.
        filter_structures (bool): Flag indicating whether to filter structure files or not.
        proteins_flag (bool): Flag - True if proteins-compare is run

    Returns:
        None
    """

    # read in amino-acid sequences
    sequences_aa = {}
    for record in SeqIO.parse(fasta_aa, "fasta"):
        sequences_aa[record.id] = str(record.seq)

    # Index structure files by CDS id (filename stem) so the per-CDS lookup
    # below is O(1). The old O(K) linear scan inside an O(K) outer loop was
    # the dominant cost on big genomes: 50k CDS × 50k files = 2.5e9 string
    # equality checks. Single-pass index keeps both ".pdb" and ".cif" hits
    # under the same key so the same len-check / take-first logic works.
    structures_by_cds_id: dict = {}
    for file in os.listdir(structure_dir):
        if file.endswith(".pdb") or file.endswith(".cif"):
            stem = file[:-4]  # ".pdb" and ".cif" are both 4 chars
            structures_by_cds_id.setdefault(stem, []).append(file)

    num_structures = 0

    # Checks that ID is in the pdbs

    no_structure_cds_ids = []

    for id in sequences_aa.keys():
        if proteins_flag:
            # will just be the CDS id if it is proteins-compare
            cds_id = id
        else:
            # Header format is "contig_id:cds_id"; everything after the first
            # colon is the CDS id (inner colons are preserved via split max=1).
            parts = id.split(":", 1)
            if len(parts) == 1:
                logger.warning(
                    f"FASTA header '{id}' has no ':' separator — expected "
                    "'contig_id:cds_id'. Treating the whole id as cds_id; "
                    "no matching structure will be found."
                )
                cds_id = id
            else:
                cds_id = parts[1]

        # record_id = id.split(":")[0]
        # this is potentially an issue if a contig has > 9999 AAs
        # need to fix with Pharokka possibly. Unlikely to occur but might!
        # enforce names as  "{cds_id}.pdb" or "{cds_id}.cif" (AF3)

        matching_files = structures_by_cds_id.get(cds_id, [])

        # delete the copying upon release, but for now do the copying to easy get the > Oct 2021 PDBs
        # with filter_structures
        if len(matching_files) == 1:
            if filter_structures is True:
                source_path = Path(structure_dir) / matching_files[0]
                destination_path = Path(filtered_structures_path) / matching_files[0]
                shutil.copyfile(source_path, destination_path)
            num_structures += 1

        # should neve happen but in case
        if len(matching_files) > 1:
            logger.warning(f"More than 1 structures found for {cds_id}")
            logger.warning("Taking the first one")
            if filter_structures is True:
                source_path = Path(structure_dir) / matching_files[0]
                destination_path = Path(filtered_structures_path) / matching_files[0]
                shutil.copyfile(source_path, destination_path)
            num_structures += 1
        elif len(matching_files) == 0:
            logger.warning(f"No structure found for {cds_id}")
            logger.warning(f"{cds_id} will be ignored in annotation")
            no_structure_cds_ids.append(cds_id)

    if num_structures == 0:
        logger.error(
            f"No structures with matching CDS ids were found at all. Check the {structure_dir} directory"
        )

    # generate the db
    short_db_name = f"{prefix}"
    structure_db_name: Path = Path(foldseek_db_path) / short_db_name
    query_structure_dir = structure_dir

    # choose the filtered directory if true
    # otherwise all pdbs in the structure_dir will be made into the foldseek DB
    if filter_structures is True:
        query_structure_dir = filtered_structures_path

    foldseek_createdb_from_structures = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createdb {query_structure_dir} {structure_db_name} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_createdb_from_structures)


#### foldseek_gpu


def create_foldseek_prostt5_gpu_db(
    fasta_aa: Path, foldseek_db_path: Path, db_dir: Path, logdir: Path
) -> None:
    """
    Convert a Foldseek DB with ProstT5 3Di predictions using Foldseek-GPU

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        foldseek_db_path (Path): Path to the directory where Foldseek database will be stored.
        db_dir (Path): Path to the Phold DB
        logdir (Path): Path to the directory where logs will be stored.
    Returns:
        None
    """

    prostt5_db_path = Path(db_dir) / "prostt5_weights"

    foldseek_createdb_prostt5 = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createdb {fasta_aa} {foldseek_db_path}  --prostt5-model {prostt5_db_path}  ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_createdb_prostt5)
