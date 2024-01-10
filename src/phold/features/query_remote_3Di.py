#!/usr/bin/env python3

from pathlib import Path
from loguru import logger
import requests
import time
import re


def query_remote_3di(cds_dict: dict, out_path: Path):
    predictions = {}

    url_base = "https://3di.foldseek.com/predict/"

    for record_id, cds_records in cds_dict.items():
        # instantiate the nested dict
        predictions[record_id] = {}

        seq_record_dict = cds_dict[record_id]
        seq_dict = {}

        # gets the seq_dict with key for id and the translation
        for cds_id, seq_feature in seq_record_dict.items():
            # get the amino acid seq
            aa_seq = seq_feature.qualifiers["translation"][0]
            seq_len = len(aa_seq)
            print(aa_seq)
            print(cds_id)

            url = f"{url_base}{aa_seq}"
            print(url)

            # Send a GET request to the URL
            response = requests.get(url)

            # Check the response status code to ensure the request was successful
            if response.status_code == 200:
                seq_3Di = response.text
                seq_3Di = str(seq_3Di)
                # it will remove any characters other than a to z, A to Z
                seq_3Di = re.sub("[^a-z]+", "", seq_3Di, flags=re.IGNORECASE)
                # Removing quotation marks
                print(seq_3Di)

            else:
                logger.warning(
                    f"Request for {cds_id} L={seq_len} failed with status code {response.status_code}"
                )

            # Add a small delay between requests
            time.sleep(0.25)

            # add the prediction to the seq
            predictions[record_id][cds_id] = seq_3Di

    # save to file
    output_3di: Path = Path(out_path) / "output3di.fasta"

    with open(output_3di, "w+") as out_f:
        for contig_id, rest in predictions.items():
            contig_predictions_dict = predictions[contig_id]

            # writes the CDS to file
            for seq_id, seq_3di in contig_predictions_dict.items():
                out_f.write(f">{contig_id}:{seq_id}\n")
                out_f.write(f"{seq_3di}\n")

    logger.info(f"Finished writing results to {out_path}")


# # Send a GET request to the URL
# response = requests.get(url)

# # Check the response status code to ensure the request was successful
# if response.status_code == 200:
#     # Parse the response content as text (JSON, HTML, etc.)
#     data = response.text
#     print(data)
# else:
#     print(f"Request failed with status code {response.status_code}")


#     output_3di: Path = Path(out_path) / "output3di.fasta"
