#!/usr/bin/env python3

import re
import time
from pathlib import Path
from typing import Dict

import requests
from loguru import logger


def query_remote_3di(
    cds_dict: Dict[str, dict], fasta_3di: Path, fasta_flag: bool
) -> None:
    """

    Query remote Foldseek ProstT5 server for 3Di predictions of amino acid sequences and write to file.

    Args:
        cds_dict (Dict[str, dict]): Dictionary containing CDS sequences.
        fasta_3di (Path): Path to save the generated 3Di FASTA file.
        fasta_flag (bool): True if input is a FASTA file

    Returns:
        None
    """

    predictions = {}

    logger.info("Querying the Foldseek ProstT5 server")
    logger.info(
        "Each CDS will be queried with a 2 second delay, so as not to overwhelm the server"
    )
    logger.info("Please be patient")

    url_base = "https://3di.foldseek.com/predict/"

    for record_id, cds_records in cds_dict.items():
        # instantiate the nested dict
        predictions[record_id] = {}

        seq_record_dict = cds_dict[record_id]
        seq_dict = {}

        # gets the seq_dict with key for id and the translation
        for cds_id, seq_feature in seq_record_dict.items():
            logger.info(f"Querying {cds_id}")
            # get the amino acid seq
            aa_seq = seq_feature.qualifiers["translation"]
            seq_len = len(aa_seq)

            url = f"{url_base}{aa_seq}"

            # Send a GET request to the URL
            response = requests.get(url)

            # Check the response status code to ensure the request was successful
            if response.status_code == 200:
                seq_3Di = response.text
                seq_3Di = str(seq_3Di)
                # it will remove any characters other than a to z, A to Z
                seq_3Di = re.sub("[^a-z]+", "", seq_3Di, flags=re.IGNORECASE)
                # Removing quotation marks

            else:
                logger.warning(
                    f"Request for {cds_id} L={seq_len} failed with status code {response.status_code}"
                )

            # Add 2s delay between requests so don't destroy the server
            # 1s delay risks https error
            time.sleep(2)

            # add the prediction to the seq
            predictions[record_id][cds_id] = seq_3Di

    # save to file
    with open(fasta_3di, "w+") as out_f:
        for contig_id, rest in predictions.items():
            contig_predictions_dict = predictions[contig_id]

            # writes the CDS to file
            for seq_id, seq_3di in contig_predictions_dict.items():
                out_f.write(f">{contig_id}:{seq_id}\n")
                out_f.write(f"{seq_3di}\n")

    logger.info(f"Finished writing results to {fasta_3di}")
