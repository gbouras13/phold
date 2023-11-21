"""
Module for manipulating genbank files
taken from phynteny
"""
import binascii
import gzip
import random
import re

# imports
import pandas as pd
from Bio import SeqIO
from loguru import logger
from pandas.errors import EmptyDataError
from pathlib import Path


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

    if is_gzip_file(genbank.strip()):
        try:
            with gzip.open(genbank.strip(), "rt") as handle:
                gb_dict = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
            handle.close()
        except ValueError:
            logger.error(f"{genbank.strip()} is not a genbank file!")
            raise

    else:
        try:
            with open(genbank.strip(), "rt") as handle:
                gb_dict = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
            handle.close()
        except ValueError:
            logger.error(f"{genbank.strip()} is not a genbank file!")
            raise

    return gb_dict


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


def extract_features(this_phage):
    """
    Extract the required features and format as a dictionary

    param this_phage: phage genome extracted from genbank file
    return: dictionary with the features for this specific phage
    """

    phage_length = len(this_phage.seq)
    this_CDS = [i for i in this_phage.features if i.type == "CDS"]  # coding sequences

    position = [
        (int(this_CDS[i].location.start), int(this_CDS[i].location.end))
        for i in range(len(this_CDS))
    ]
    sense = [
        re.split("]", str(this_CDS[i].location))[1][1] for i in range(len(this_CDS))
    ]
    protein_id = [
        this_CDS[i].qualifiers.get("protein_id") for i in range(len(this_CDS))
    ]
    protein_id = [p[0] if p is not None else None for p in protein_id]
    phrogs = [this_CDS[i].qualifiers.get("phrog") for i in range(len(this_CDS))]
    phrogs = ["No_PHROG" if i is None else i[0] for i in phrogs]

    return {
        "length": phage_length,
        "phrogs": phrogs,
        "protein_id": protein_id,
        "sense": sense,
        "position": position,
    }


def filter_mmseqs(phrog_output, Eval=1e-5):
    """
    Function to filter the phogs mmseqs output
    If there are two equally good annotations, the annotation with the greatest coverage is taken, then e

    param phrog_output: dataframe of phrog annotations
    param Eval: evalue to filter annotations default (1e-5)
    return: dictionary of the phrog annotations
    """

    # rename the headers
    phrog_output.columns = [
        "phrog",
        "seq",
        "alnScore",
        "seqIdentity",
        "eVal",
        "qStart",
        "qEnd",
        "qLen",
        "tStart",
        "tEnd",
        "tLen",
    ]
    phrog_output["coverage"] = phrog_output["tEnd"] - phrog_output["tStart"] + 1
    phrog_output["phrog"] = [
        re.split("_", p)[1] for p in phrog_output["phrog"]
    ]  # convert phrog to a number
    phrog_output["phrog"] = [re.split("#", p)[0][:-1] for p in phrog_output["phrog"]]

    # filter to have an e-value lower than e-5 - can change this not to be hardcoded
    phrog_output = phrog_output[phrog_output["eVal"] < Eval]

    # filter annotations with multiple hits
    phrog_output = (
        phrog_output.groupby("seq", as_index=False).coverage.max().merge(phrog_output)
    )  # hit with greatest coverage
    phrog_output = (
        phrog_output.groupby("seq", as_index=False).eVal.min().merge(phrog_output)
    )  # hit with the best evalue
    phrog_output = (
        phrog_output.groupby("seq", as_index=False).qLen.min().merge(phrog_output)
    )  # hit with the shortest query length

    return dict(zip(phrog_output["seq"].values, phrog_output["phrog"].values))


def add_predictions(gb_dict, predictions):
    """
    Add predictions to the genbank dictionary

    param gb_dict: genbank file as a dictionary
    param predictions: predictions to add to the genbank file
    return updated dictionary with features
    """

    keys = list(gb_dict.keys())

    for i in range(len(predictions)):
        gb_dict[keys[i]]["phynteny"] = predictions[i]
    return gb_dict


def write_genbank(gb_dict, filename):
    """
    write genbank dictionary to a file
    """

    keys = list(gb_dict.keys())

    # check for gzip
    if filename.strip()[-3:] == ".gz":
        with gzip.open(filename, "wt") as handle:
            for key in keys:
                SeqIO.write(gb_dict.get(key), handle, "genbank")
        handle.close()

    else:
        with open(filename, "wt") as handle:
            for key in keys:
                SeqIO.write(gb_dict.get(key), handle, "genbank")
        handle.close()
