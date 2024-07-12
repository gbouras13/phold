#!/usr/bin/env python3
from pathlib import Path

import h5py
import pandas as pd
from Bio import SeqIO


def split_3Di_embeddings_by_prob(
    fasta_aa: Path,
    fasta_3di: Path,
    embeddings_h5: Path,
    low_prob_embeddings_h5: Path,
    probs_3di: Path,
    output: Path,
    split_threshold: float,
) -> None:
    """
    Split embeddings based on 3Di probability thresholds from ProstT5.

    Args:
        fasta_aa (Path): Path to the amino acid FASTA file.
        fasta_3di (Path): Path to the 3DI FASTA file.
        embeddings_h5 (Path): Path to the per protein embeddings h5 file.
        low_prob_embeddings_h5 (Path): Path to the per protein low probaility embeddings h5 file.
        probs_3di (Path): Path to the CSV file containing probabilities associated with 3DI sequences.
        output (Path): Path to the directory where output files will be saved.
        split_threshold (float): Probability threshold for splitting sequences.

    Returns:
        None
    """

    probs_3di_df = pd.read_csv(
        probs_3di, header=None, names=["cds_id", "prostt5_prob"], sep=","
    )

    # Create two sets based on the split_threshold
    high_prob_set = set(
        probs_3di_df[probs_3di_df["prostt5_prob"] >= split_threshold]["cds_id"]
    )
    low_prob_set = set(
        probs_3di_df[probs_3di_df["prostt5_prob"] < split_threshold]["cds_id"]
    )

    # write the high pro 3dis out
    high_prob_fasta_3di_out_path: Path = Path(output) / "high_prob_3di.fasta"

    # Open output files for writing
    high_prob_3di_out_file = open(high_prob_fasta_3di_out_path, "w")

    # Open and read the FASTA file
    with open(fasta_3di, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract cds_id from the header
            cds_id = record.id.split(":")[1]

            # Write the record to the appropriate output file based on the cds_id
            if cds_id in high_prob_set:
                SeqIO.write(record, high_prob_3di_out_file, "fasta")

    # Close the output files
    high_prob_3di_out_file.close()

    # write the aas out
    high_prob_fasta_aa_out_path: Path = Path(output) / "high_prob_aa.fasta"
    low_prob_fasta_aa_out_path: Path = Path(output) / "low_prob_aa.fasta"

    # Open output files for writing
    high_prob_aa_out_file = open(high_prob_fasta_aa_out_path, "w")
    low_prob_aa_out_file = open(low_prob_fasta_aa_out_path, "w")

    # Open and read the FASTA file
    with open(fasta_aa, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract cds_id from the header
            cds_id = record.id.split(":")[1]

            # Write the record to the appropriate output file based on the cds_id
            if cds_id in high_prob_set:
                SeqIO.write(record, high_prob_aa_out_file, "fasta")
            else:
                SeqIO.write(record, low_prob_aa_out_file, "fasta")

    # Close the output files
    high_prob_aa_out_file.close()
    low_prob_aa_out_file.close()

    # embeddings now

    # Open the HDF5 file for reading
    with h5py.File(low_prob_embeddings_h5, "w") as hfw:
        with h5py.File(embeddings_h5, "r") as hfr:

            for embeddings_id in hfr.keys():
                # to get cds id - need the contig to be saved in the embeddings (for downstream processing ease)
                cds_id = embeddings_id.split(":")[1]

                if cds_id in low_prob_set:
                    dataset = hfr[embeddings_id]
                    hfw.create_dataset(embeddings_id, data=dataset)


# def split_3di_fasta_by_prob(
#     fasta_aa: Path,
#     fasta_3di: Path,
#     probs_3di: Path,
#     output: Path,
#     split_threshold: float,
# ) -> None:
#     """
#     Split sequences in 3Di FASTA files based on probability thresholds.

#     Args:
#         fasta_aa (Path): Path to the amino acid FASTA file.
#         fasta_3di (Path): Path to the 3DI FASTA file.
#         probs_3di (Path): Path to the CSV file containing probabilities associated with 3DI sequences.
#         output (Path): Path to the directory where output files will be saved.
#         split_threshold (float): Probability threshold for splitting sequences.

#     Returns:
#         None
#     """


#     probs_3di_df = pd.read_csv(
#         probs_3di, header=None, names=["cds_id", "prostt5_prob"], sep=","
#     )

#     # Create two sets based on the split_threshold
#     high_prob_set = set(
#         probs_3di_df[probs_3di_df["prostt5_prob"] >= split_threshold]["cds_id"]
#     )
#     low_prob_set = set(
#         probs_3di_df[probs_3di_df["prostt5_prob"] < split_threshold]["cds_id"]
#     )

#     # write the 3dis out

#     high_prob_fasta_3di_out_path: Path = Path(output) / "high_prob_3di.fasta"
#     low_prob_fasta_3di_out_path: Path = Path(output) / "low_prob_3di.fasta"

#     # Open output files for writing
#     high_prob_3di_out_file = open(high_prob_fasta_3di_out_path, "w")
#     low_prob_3di_out_file = open(low_prob_fasta_3di_out_path, "w")

#     # Open and read the FASTA file
#     with open(fasta_3di, "r") as fasta_file:
#         for record in SeqIO.parse(fasta_file, "fasta"):
#             # Extract cds_id from the header
#             cds_id = record.id.split(":")[1]

#             # Write the record to the appropriate output file based on the cds_id
#             if cds_id in high_prob_set:
#                 SeqIO.write(record, high_prob_3di_out_file, "fasta")
#             elif cds_id in low_prob_set:
#                 SeqIO.write(record, low_prob_3di_out_file, "fasta")

#     # Close the output files
#     high_prob_3di_out_file.close()
#     low_prob_3di_out_file.close()

#     # write the aas out

#     high_prob_fasta_aa_out_path: Path = Path(output) / "high_prob_aa.fasta"
#     low_prob_fasta_aa_out_path: Path = Path(output) / "low_prob_aa.fasta"

#     # Open output files for writing
#     high_prob_aa_out_file = open(high_prob_fasta_aa_out_path, "w")
#     low_prob_aa_out_file = open(low_prob_fasta_aa_out_path, "w")

#     # Open and read the FASTA file
#     with open(fasta_aa, "r") as fasta_file:
#         for record in SeqIO.parse(fasta_file, "fasta"):
#             # Extract cds_id from the header
#             cds_id = record.id.split(":")[1]

#             # Write the record to the appropriate output file based on the cds_id
#             if cds_id in high_prob_set:
#                 SeqIO.write(record, high_prob_aa_out_file, "fasta")
#             elif cds_id in low_prob_set:
#                 SeqIO.write(record, low_prob_aa_out_file, "fasta")

#     # Close the output files
#     high_prob_aa_out_file.close()
#     low_prob_aa_out_file.close()
