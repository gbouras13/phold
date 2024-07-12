#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@authors: gbouras13 and @mheinzinger 

Code adapted from @mheinzinger 

https://github.com/Rostlab/EAT

This provides class and functions to use ProstT5 embeddings to do EAT for low confidence ProstT5 CNN predictions

"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import h5py
import numpy as np
import torch
import torch.nn.functional as F
from loguru import logger


# EAT: Embedding-based Annotation Transfer
class EAT:
    """
    Taken from https://github.com/Rostlab/EAT/eat.py with modifications
    """

    def __init__(
        self, lookup_path, query_path, output_dir, num_NN, cpu, split_threshold
    ):

        global device

        if cpu is True:
            device = torch.device("cpu")
            dev_name = "cpu"
        else:
            # check for NVIDIA/cuda
            if torch.cuda.is_available():
                device = torch.device("cuda:0")
                dev_name = "cuda:0"
            # check for apple silicon/metal
            elif torch.backends.mps.is_available():
                device = torch.device("mps")
                dev_name = "mps"
            else:
                device = torch.device("cpu")
                dev_name = "cpu"

        # logger device only if the function is called
        logger.info("Using device: {} for EAT".format(dev_name))

        self.num_NN = num_NN
        self.Embedder = None

        self.output_dir = output_dir

        self.query_ids, self.query_embs = self.read_inputs(query_path, split_threshold)
        logger.info("Query embeddings read")
        self.lookup_ids, self.lookup_embs = self.read_inputs(
            lookup_path, split_threshold
        )
        logger.info("Phold DB embeddings read")

    def read_inputs(self, lookup_path, split_threshold):
        """
        has to be a .h5 file
        """

        # define path for storing embeddings
        if not lookup_path.is_file():
            logger.warning("No embedding H5 could be found for: {}".format(lookup_path))
            logger.error("Files are expected to either end with .fasta or .h5. Exiting")
            raise FileNotFoundError

        if lookup_path.name.endswith(".h5"):  # if the embedding file already exists
            return self.read_embeddings(lookup_path, split_threshold)
        else:
            logger.error(
                "The file you passed neither ended with .fasta nor .h5. Only those file formats are currently supported."
            )
            raise NotImplementedError

    def read_embeddings(self, emb_path, split_threshold):

        # read in the embeddings
        h5_f = h5py.File(emb_path, "r")
        # get dataset

        dataset = {pdb_id: np.array(embd) for pdb_id, embd in h5_f.items()}

        if len(dataset.items()) == 0:
            logger.error(
                f"No proteins had ProstT5 confidence under --split_threshold {split_threshold}. Please check this or do not use --split."
            )
        keys, embeddings = zip(*dataset.items())
        # if keys[0].startswith("cath"):
        #     keys = [key.split("|")[2].split("_")[0] for key in keys ]
        # matrix of values (protein-embeddings); n_proteins x embedding_dim
        # get embeddings
        embeddings = np.vstack(embeddings)
        return list(keys), torch.tensor(embeddings).to(device).float()

    # # Open the HDF5 file for reading
    # with h5py.File(str(out_path), "r") as hf:
    #     # Iterate over each group (contig_id)
    #     for contig_id in hf.keys():
    #         # Get the group corresponding to the contig_id
    #         contig_group = hf[contig_id]

    #         # Iterate over each dataset (sequence_id) in the group
    #         for sequence_id in contig_group.keys():
    #             # Get the dataset corresponding to the sequence_id
    #             dataset = contig_group[sequence_id]

    #             # Read the data from the dataset
    #             embedding_data = dataset[:]

    #             # Now you have the embedding data for this sequence_id
    #             # You can process it as needed
    #             print(f"Contig ID: {contig_id}, Sequence ID: {sequence_id}, Embedding: {embedding_data}")

    def pdist(self, lookup, queries, norm=2, use_double=False):
        lookup = lookup.unsqueeze(dim=0)
        queries = queries.unsqueeze(dim=0)
        # double precision improves performance slightly but can be removed for speedy predictions (no significant difference in performance)
        if use_double:
            lookup = lookup.double()
            queries = queries.double()

        torch.cuda.empty_cache()

        try:  # try to batch-compute all pairwise-distance on GPU
            pdist = torch.cdist(lookup, queries, p=norm).squeeze(dim=0)
        except RuntimeError as e:
            logger.warning("Encountered RuntimeError: {}".format(e))
            # next try computing in batches of batch_size
            try:
                batch_size = 200
                logger.warning(
                    f"Trying batched inference on GPU with batchsize {batch_size}."
                )

                results = []

                for start_idx in range(0, queries.shape[1], batch_size):
                    end_idx = min(start_idx + batch_size, queries.shape[1])

                    with torch.no_grad():
                        # Process a batch of queries
                        query_batch = queries[0:1, start_idx:end_idx]
                        # Compute the distance for the current batch
                        dist_batch = torch.cdist(lookup, query_batch, p=norm).squeeze(
                            dim=0
                        )
                        results.append(
                            dist_batch.cpu()
                        )  # Move to CPU to free up GPU memory

                    # Clear CUDA cache to free up memory
                    torch.cuda.empty_cache()

                # Stack all the results
                pdist = torch.cat(results, dim=1)

            except RuntimeError as e:
                logger.warning("Encountered RuntimeError: {}".format(e))
                logger.warning(
                    "Trying single query inference on GPU. This may take a while."
                )
                try:
                    # michael's original code -> causes memory leakage
                    # therefore need to loop over and free GPU mem
                    # if OOM for batch-GPU, re-try single query pdist computation on GPU
                    # pdist = (
                    #     torch.stack(
                    #         [
                    #             torch.cdist(lookup, queries[0:1, q_idx], p=norm).squeeze(
                    #                 dim=0
                    #             )
                    #             for q_idx in range(queries.shape[1])
                    #         ]
                    #     )
                    #     .squeeze(dim=-1)
                    #     .T
                    # )

                    # to avoid OOM -> memory build-up

                    results = []

                    for q_idx in range(queries.shape[1]):
                        with torch.no_grad():
                            # Compute the distance for the current query
                            dist = torch.cdist(
                                lookup, queries[0:1, q_idx], p=norm
                            ).squeeze(dim=0)
                            # Append the result to the list
                            results.append(
                                dist.cpu()
                            )  # Move to CPU to free up GPU memory

                        # Clear CUDA cache to free up memory
                        torch.cuda.empty_cache()

                    # Stack the results and perform necessary squeezing and transposing
                    pdist = torch.stack(results).squeeze(dim=-1).T

                except (
                    RuntimeError
                ) as e:  # if OOM for single GPU, re-try single query on CPU
                    logger.warning("Encountered RuntimeError: {}".format(e))
                    logger.warning("Trying to move single query computation to CPU.")
                    lookup = lookup.to("cpu")
                    queries = queries.to("cpu")
                    pdist = (
                        torch.stack(
                            [
                                torch.cdist(
                                    lookup, queries[0:1, q_idx], p=norm
                                ).squeeze(dim=0)
                                for q_idx in range(queries.shape[1])
                            ]
                        )
                        .squeeze(dim=-1)
                        .T
                    )

        return pdist

    def get_NNs(self, threshold, random=False):
        p_dist = self.pdist(self.lookup_embs, self.query_embs)

        nn_dists, nn_idxs = torch.topk(p_dist, self.num_NN, largest=False, dim=0)

        nn_dists, nn_idxs = nn_dists.to("cpu"), nn_idxs.to("cpu")
        predictions = list()
        n_test = len(self.query_ids)
        for test_idx in range(n_test):  # for all test proteins
            query_id = self.query_ids[test_idx]  # get id of test protein
            nn_idx = nn_idxs[:, test_idx]
            nn_dist = nn_dists[:, test_idx]
            for nn_iter, (nn_i, nn_d) in enumerate(zip(nn_idx, nn_dist)):
                # index of nearest neighbour (nn) in train set
                nn_i, nn_d = int(nn_i), float(nn_d)

                # if a threshold is passed, skip all proteins above this threshold
                if threshold is not None and nn_d > threshold:
                    continue
                # get id of nn (infer annotation)
                lookup_id = self.lookup_ids[nn_i]
                # lookup_label = self.lookupLabels[lookup_id]
                # query_label = self.queryLabels[query_id]
                predictions.append((query_id, lookup_id, nn_d, nn_iter))
        logger.info("Computing embedding nearest neighbours finished")
        return predictions

    def write_predictions(self, predictions, prefix):
        out_p = self.output_dir / f"{prefix}_eat_result.txt"
        with open(out_p, "w+") as out_f:
            out_f.write(
                "Query-ID\tLookup-ID\tEmbedding distance\tNearest-Neighbor-Idx\n"
            )
            out_f.write(
                "\n".join(
                    [
                        "{}\t{}\t{:.4f}\t{}".format(
                            query_id, lookup_id, eat_dist, nn_iter + 1
                        )
                        for query_id, lookup_id, eat_dist, nn_iter in predictions
                    ]
                )
            )

        # dummy temp tsv for foldseek
        out_p = self.output_dir / "foldseek_results_low.tsv"

        with open(out_p, "w+") as out_f:
            for query_id, lookup_id, eat_dist, nn_iter in predictions:
                out_f.write(f"{query_id}\t{lookup_id}\t0\t0\t0\t0\t0\t0\t0\t0\t0\n")

        # col_list = [
        #     "query",
        #     "target",
        #     "bitscore",
        #     "fident",
        #     "evalue",
        #     "qStart",
        #     "qEnd",
        #     "qLen",
        #     "tStart",
        #     "tEnd",
        #     "tLen",
        # ]

        return None


def run_eat(
    out_path: Path,
    prefix: str,
    db_dir: Path,
    output_h5_per_protein: Path,
    cpu: bool = False,
    proteins_flag: bool = False,
    num_NN: int = 1,
    eat_threshold: float = 1.00,
    split_threshold: float = 60,
) -> bool:
    """
    Run EAT against on embeddings using euclidean distance

    Args:
        out_path (Path): Path to the output directory.
        prefix (str): Prefix for the output files.
        db_dir (Path): Directory containing the pre-calculated ProstT5 embeddings.
        output_h5_per_protein (Path): Path to the output h5 per proteins embeddings file.
        proteins_flag (bool, optional): Whether the sequences are proteins. Defaults to False.
        num_NN (int): Number of nearest neighbours to search for. Defaults to 1
        eat_threshold (float): Threshold of euclidean distance to search for EAT.
        split_threshold (float): Threshold for splitting EAT vs Foldseek --split_threshold


    Returns:
        bool: True if embeddings and predictions are generated successfully.
    """

    # db
    lookup_p = Path(db_dir) / "phold_db.h5"
    # query
    query_p = Path(output_h5_per_protein)

    # output directory
    output_d = Path(out_path)

    assert num_NN > 0, logger.error(
        "Only positive number of nearest neighbors can be retrieved."
    )

    eater = EAT(lookup_p, query_p, output_d, num_NN, cpu, split_threshold)
    logger.info("Computing embeddings nearest neighbours")
    predictions = eater.get_NNs(threshold=eat_threshold)
    eater.write_predictions(predictions, prefix)

    return None
