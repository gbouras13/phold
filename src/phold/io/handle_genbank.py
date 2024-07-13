"""
Module for manipulating genbank files
some taken from phynteny https://github.com/susiegriggo/Phynteny
"""

import binascii
import gzip
import multiprocessing.pool
from datetime import datetime
from pathlib import Path
from typing import Dict

import pandas as pd
import pyrodigal_gv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

# imports


def is_gzip_file(f: Path) -> bool:
    """
    Method copied from Phispy see https://github.com/linsalrob/PhiSpy/blob/master/PhiSpyModules/helper_functions.py

    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    Args:
        f (Path): The file to test.

    Returns:
        bool: True if the file is gzip compressed, otherwise False.
    """
    with open(f, "rb") as i:
        return binascii.hexlify(i.read(2)) == b"1f8b"


def get_genbank(genbank: Path) -> dict:
    """
    Convert a GenBank file to a dictionary.

    This function reads a GenBank file and converts it into a dictionary.

    Args:
        genbank (Path): Path to the GenBank file.

    Returns:
        dict: A dictionary representation of the GenBank file.

    Raises:
        ValueError: If the provided file is not a GenBank file.
    """

    logger.info(f"Checking if input {genbank} is a Genbank file")

    if is_gzip_file(genbank.strip()):
        try:
            with gzip.open(genbank.strip(), "rt") as handle:
                gb_dict = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
            handle.close()
        except ValueError:
            logger.warning(f"{genbank.strip()} is not a genbank file")
            raise

    else:
        try:
            with open(genbank.strip(), "rt") as handle:
                gb_dict = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
            handle.close()
        except ValueError:
            logger.warning(f"{genbank} is not a genbank file")
            raise

    return gb_dict


def get_fasta_run_pyrodigal_gv(input: Path, threads: int) -> dict:
    """

    Check if a file is in the nucleotide FASTA format. If so, run pyrodigal-gv and convert the CDS to a dictionary.

    Args:
        input (Path): Path to the FASTA file.

    Returns:
        dict: A dictionary representation of the CDS in the FASTA file.

    Raises:
        ValueError: If the provided file is not a FASTA file.
    """

    if is_gzip_file(input.strip()):
        try:
            with gzip.open(input.strip(), "rt") as handle:
                # gb_dict = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        except ValueError:
            logger.warning(f"{input} is not a FASTA file")
            logger.error(
                f"Your input {input} is neither Genbank nor FASTA format. Please check your input"
            )
            raise
    else:
        try:
            with open(input.strip(), "rt") as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        except ValueError:
            logger.warning(f"{input} is not a FASTA file")
            logger.error(
                f"Your input {input} is neither Genbank nor FASTA format. Please check your input"
            )
            raise

    # then run pyrodigal

    gb_dict = {}

    orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)

    def _find_genes(record):
        genes = orf_finder.find_genes(str(record.seq))
        return (record.id, record.seq, genes)

    def run_pool(pool, records):
        for record_id, record_seq, genes in pool.imap(_find_genes, records):
            i = 0
            all_features = []
            for gene in genes:
                i += 1
                location = FeatureLocation(
                    start=gene.begin, end=gene.end, strand=gene.strand
                )
                feature = SeqFeature(location, type="CDS")
                counter = "{:04d}".format(i)
                cds_id = f"{record_id}_CDS_" + counter
                feature.qualifiers["ID"] = cds_id
                feature.qualifiers["function"] = "unknown function"
                feature.qualifiers["product"] = "hypothetical protein"
                feature.qualifiers["phrog"] = "No_PHROG"
                feature.qualifiers["source"] = (
                    f"Pyrodigal-gv_{pyrodigal_gv.__version__}"
                )
                feature.qualifiers["transl_table"] = gene.translation_table
                # from the API
                # translation_table (int, optional) – An alternative translation table to use to translate the gene.
                # Use None (the default) to translate using the translation table this gene was found with.
                feature.qualifiers["translation"] = gene.translate(
                    include_stop=False
                ).upper()
                all_features.append(feature)

            seq_record = SeqIO.SeqRecord(
                seq=Seq(record_seq), id=record_id, description="", features=all_features
            )
            gb_dict[record_id] = seq_record

        return gb_dict

    with multiprocessing.pool.ThreadPool(threads) as pool:
        if is_gzip_file(input.strip()):
            with gzip.open(input.strip(), "rt") as handle:
                records = SeqIO.parse(handle, "fasta")
                gb_dict = run_pool(pool, records)
        else:
            with open(input.strip(), "rt") as handle:
                records = SeqIO.parse(handle, "fasta")
                gb_dict = run_pool(pool, records)

    return gb_dict


def write_genbank(
    updated_cds_dict: Dict[str, Dict[str, SeqFeature]],
    non_cds_dict: Dict[str, Dict[str, SeqFeature]],
    source_dict: Dict[str, Dict[str, SeqFeature]],
    prefix: str,
    gb_dict: Dict[str, SeqRecord],
    output: Path,
    proteins_flag: bool,
    separate: bool,
    fasta_flag: bool,
) -> pd.DataFrame:
    """
    Write sequences to GenBank files.

    Args:
        updated_cds_dict (Dict[str, Dict[str, SeqFeature]]): Dictionary containing updated CDS features.
        non_cds_dict (Dict[str, Dict[str, SeqFeature]]): Dictionary containing non-CDS features.
        source_dict (Dict[str, Dict[str, SeqFeature]]): Dictionary containing CDS source assignments.
        prefix (str): Prefix for the output GenBank file.
        gb_dict (Dict[str, SeqRecord]): Dictionary containing SeqRecords.
        output (Union[str, Path]): Path to the output directory.
        proteins_flag (bool): Flag indicating whether proteins are present.
        separate (bool): Flag indicating whether to write separate GenBank files.
        fasta_flag (bool): Flag indicating whether input is a FASTA file.

    Returns:
        pd.DataFrame: DataFrame containing information about each CDS.
    """

    # separate gbks per contig
    if separate is True:
        separate_output = Path(output) / "separate_gbks"
        separate_output.mkdir(parents=True, exist_ok=True)

    seq_records = []
    per_cds_list = []

    for record_id, record in gb_dict.items():
        # Merge updated_cds_dict and non_cds_dict
        merged_dict = {
            record_id: {
                **updated_cds_dict.get(record_id, {}),
                **non_cds_dict.get(record_id, {}),
            }
        }

        # Extract features into a list
        all_features = [
            feature
            for features in merged_dict.values()
            for feature in features.values()
        ]

        # Sort features based on the beginning position
        sorted_features = sorted(all_features, key=lambda x: x.location.start)

        # clean cds_feature and append for dataframe
        for cds_feature in sorted_features:
            if cds_feature.type == "CDS":
                if proteins_flag is True:
                    cds_info = {
                        "cds_id": cds_feature.qualifiers["ID"],
                        "phrog": cds_feature.qualifiers["phrog"][0],
                        "function": cds_feature.qualifiers["function"][0],
                        "product": cds_feature.qualifiers["product"][0],
                    }
                else:
                    # because for some reason when parsing the pharokka genbank, it is a list, fasta it is not
                    if fasta_flag is True:
                        try:
                            transl_table = cds_feature.qualifiers["transl_table"]
                        except:
                            # for older pharokka input before v1.5.0
                            transl_table = "11"
                    else:
                        try:
                            transl_table = cds_feature.qualifiers["transl_table"][0]
                        except:
                            # for older pharokka input before v1.5.0
                            transl_table = "11"

                    # to reverse the start and end coordinates for output tsv + fix genbank 0 index start relative to pharokka
                    if cds_feature.location.strand == -1:  # neg strand
                        start = cds_feature.location.end
                        if fasta_flag is True:  # pyrodigal
                            end = cds_feature.location.start
                        else:  # genbank
                            end = cds_feature.location.start + 1
                    else:  # pos strand
                        if fasta_flag is True:  # pyrodigal
                            start = cds_feature.location.start
                        else:  # genbank
                            start = cds_feature.location.start + 1
                        end = cds_feature.location.end

                    if fasta_flag is True:
                        cds_id = cds_feature.qualifiers["ID"]
                    else:  # because for some reason when parsing the pharokka genbank, it is a list
                        cds_id = cds_feature.qualifiers["ID"][0]

                    cds_info = {
                        "contig_id": record_id,
                        "cds_id": cds_id,
                        "start": start,
                        "end": end,
                        "strand": cds_feature.location.strand,
                        "phrog": cds_feature.qualifiers["phrog"][0],
                        "function": cds_feature.qualifiers["function"][0],
                        "product": cds_feature.qualifiers["product"][0],
                        "annotation_method": source_dict[record_id][cds_id],
                        "transl_table": transl_table,
                    }

                    # Remove unwanted gbk attributes if they exist
                    keys_to_remove = [
                        "top_hit",
                        "score",
                        "phase",
                        "CARD_short_name",
                        "AMR_Gene_Family",
                        "CARD_species",
                        "vfdb_short_name",
                        "vfdb_description",
                        "vfdb_species",
                    ]
                    for key in keys_to_remove:
                        # will remove the keys
                        deleted_value = cds_feature.qualifiers.pop(key, None)
                    # get dataframe
                per_cds_list.append(cds_info)

        # write out the record to GBK file
        if proteins_flag is False:
            sequence = record.seq
        # if proteins dummy seq
        else:
            sequence = Seq("")
        seq_record = SeqIO.SeqRecord(
            seq=sequence, id=record_id, description="", features=sorted_features
        )
        seq_records.append(seq_record)

        # update the molecule type, data file division and date
        seq_record.annotations["molecule_type"] = "DNA"
        seq_record.annotations["data_file_division"] = "PHG"
        seq_record.annotations["date"] = str(
            datetime.now().strftime("%d-%b-%Y").upper()
        )

        if separate is True and proteins_flag is False:
            output_gbk_path: Path = Path(separate_output) / f"{record_id}.gbk"
            with open(output_gbk_path, "w") as output_file:
                SeqIO.write(seq_record, output_file, "genbank")

    per_cds_df = pd.DataFrame(per_cds_list)

    if proteins_flag is False:
        # convert strand
        per_cds_df["strand"] = per_cds_df["strand"].apply(
            lambda x: "-" if x == -1 else ("+" if x == 1 else x)
        )

        # only write the gbk if proteins_flag is False
        output_gbk_path: Path = Path(output) / f"{prefix}.gbk"
        with open(output_gbk_path, "w") as output_file:
            SeqIO.write(seq_records, output_file, "genbank")

    return per_cds_df


def get_proteins(fasta: Path) -> dict:
    """
    Convert an Amino Acid FASTA file to a dictionary.

    This function reads a AA FASTA file and converts it into a dictionary.

    Args:
        fasta (Path): Path to the FASTA file.

    Returns:
        dict: A dictionary representation of the FASTA file.

    Raises:
        ValueError: If the provided file is not a FASTA file.
    """

    if is_gzip_file(fasta.strip()):
        try:
            fasta_dict = {}
            with gzip.open(fasta.strip(), "rt") as handle:
                sequence_id = ""
                sequence = ""
                for line in handle:
                    line = line.strip()
                    if line.startswith(">"):
                        if sequence_id:
                            fasta_dict[sequence_id] = sequence
                        sequence_id = line[1:]
                        sequence = ""
                    else:
                        sequence += line
                if sequence_id:
                    fasta_dict[sequence_id] = sequence
            handle.close()
        except ValueError:
            logger.error(f"{fasta.strip()} is not a FASTA file!")
            raise

    else:
        try:
            fasta_dict = {}
            with open(fasta.strip(), "rt", errors="ignore") as handle:
                sequence_id = ""
                sequence = ""
                for line in handle:
                    line = line.strip()
                    if line.startswith(">"):
                        if sequence_id:
                            fasta_dict[sequence_id] = sequence
                        sequence_id = line[1:]
                        sequence = ""
                    else:
                        sequence += line
                if sequence_id:
                    fasta_dict[sequence_id] = sequence
            handle.close()
        except ValueError:
            logger.error(f"{fasta.strip()} is not a FASTA file!")
            raise

    return fasta_dict
