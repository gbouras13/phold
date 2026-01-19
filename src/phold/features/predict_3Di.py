#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Code adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/predict_3Di_encoderOnly.py

"""

import math
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from collections import defaultdict


import contextlib
import os
import h5py
import numpy as np
import torch
import torch.nn.functional as F
from loguru import logger
from torch import nn
from tqdm import tqdm
from transformers import T5EncoderModel, T5Tokenizer, AutoModel, AutoTokenizer

from phold.databases.db import check_model_download, download_zenodo_model
from phold.utils.constants import CNN_DIR



# Convolutional neural network (two convolutional layers)
class CNN(nn.Module):
    def __init__(self):
        """
        Initialize the Convolutional Neural Network (CNN) model.
        """
        super(CNN, self).__init__()

        self.classifier = nn.Sequential(
            nn.Conv2d(1024, 32, kernel_size=(7, 1), padding=(3, 0)),  # 7x32
            nn.ReLU(),
            nn.Dropout(0.0),
            nn.Conv2d(32, 20, kernel_size=(7, 1), padding=(3, 0)),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Perform forward pass through the CNN.

        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, sequence_length, embedding_size).

        Returns:
            torch.Tensor: Output tensor of shape (batch_size, sequence_length, num_classes).

        L = protein length
        B = batch-size
        F = number of features (1024 for embeddings)
        N = number of classes (20 for 3Di)
        """

        # Permute input tensor to match expected shape
        # Input shape: (batch_size, sequence_length, embedding_size)
        # Output shape: (batch_size, embedding_size, sequence_length, 1)

        x = x.permute(0, 2, 1).unsqueeze(
            dim=-1
        )  # IN: X = (B x L x F); OUT: (B x F x L, 1)

        # Pass the input through the classifier
        # Output shape: (batch_size, num_classes, sequence_length, 1)
        Yhat = self.classifier(x)  # OUT: Yhat_consurf = (B x N x L x 1)

        # Remove the singleton dimension from the output tensor
        # Output shape: (batch_size, num_classes, sequence_length)
        Yhat = Yhat.squeeze(dim=-1)  # IN: (B x N x L x 1); OUT: ( B x L x N )
        return Yhat


def get_model(
    model_dir: Path, model_name: str, cpu: bool, threads: int
) -> (AutoModel, AutoTokenizer):
    """
    Loads a HF model and tokenizer.

    Args:
        model_dir (Path): Directory where the model and tokenizer is be stored.
        model_name (str): Name of the pre-trained T5 model.
        cpu (bool): Whether to use CPU only.
        threads (int): Number of cpu threads.

    Returns:
        Tuple[AutoModel, AutoTokenizer]: Tuple containing the loaded model and tokenizer.
    """

    # sets the device

    # Torch load will map back to device from state, which often is GPU:0.
    # to overcome, need to explicitly map to active device

    global device

    torch.set_num_threads(threads)

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
            if cpu is not True:
                logger.warning(
                    "No available GPU was found, but --cpu was not specified"
                )
                logger.warning("ProstT5 will be run with CPU only")

    # logger device only if the function is called
    logger.info("Using device: {}".format(dev_name))

    # make dir if doesnt exist
    Path(model_dir).mkdir(parents=True, exist_ok=True)

    # load
    logger.info(f"Loading {model_name} from: {model_dir}/{model_name}")
    logger.info(f"If {model_dir}/{model_name} is not found, it will be downloaded")

    # check ProstT5 is downloaded
    # flag assumes transformers takes from local file (see #44)
    localfile = True
    download = False

    download = check_model_download(model_dir, model_name)

    if download is True:
        localfile = False
        logger.info(f"{model_name} not found. Downloading {model_name} from Hugging Face")
    try:
        if model_name != "gbouras13/modernprost-base":
            model = T5EncoderModel.from_pretrained(
                model_name,
                cache_dir=f"{model_dir}/",
                force_download=False,
                local_files_only=True,
            ).to(device)
        else:
            # still doesnt fully work actually, only removes some of the prints
            with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                model = AutoModel.from_pretrained(
                    model_name,
                    trust_remote_code=True,
                    cache_dir=f"{model_dir}/",
                    force_download=download,
                    local_files_only=localfile,
                ).to(device)

    except:
        logger.warning("Download from Hugging Face failed. Trying backup from Zenodo.")
        logdir = f"{model_dir}/logdir"
        if model_name != "gbouras13/modernprost-base":
            download_zenodo_model(model_dir, logdir, threads, model="ProstT5")
            model = T5EncoderModel.from_pretrained(
                model_name,
                cache_dir=f"{model_dir}/",
                force_download=False,
                local_files_only=True,
            ).to(device)
        else:
            download_zenodo_model(model_dir, logdir, threads, model="modernprost")
            with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                model = AutoModel.from_pretrained(
                model_name,
                trust_remote_code=True,
                cache_dir=f"{model_dir}/",
                force_download=download,
                local_files_only=localfile,
            ).to(device)
            model = model.half()

    model = model.eval()

    # tokeniser

    if model_name != "gbouras13/modernprost-base":
        vocab = T5Tokenizer.from_pretrained(
            model_name, cache_dir=f"{model_dir}/", do_lower_case=False
        )
    else:
        vocab = AutoTokenizer.from_pretrained(
            model_name,
            trust_remote_code=True,
            cache_dir=f"{model_dir}/",
            force_download=download,
            local_files_only=localfile,
        )


    logger.info(f"{model_name} loaded")

    return model, vocab


def write_embeddings(
    embeddings: Dict[str, Dict[str, Tuple[List[str], Any, Any]]],
    out_path: Path,
) -> None:
    """
    Write embeddings to an output file.

    Args:
        embeddings (Dict[str, Dict[str, Tuple[List[str], Any, Any]]]): Predictions dictionary containing contig IDs, sequence IDs, predictions, and additional information.
        out_path (Path): Path to the output file.

    Returns:
        None
    """

    with h5py.File(str(out_path), "w") as hf:
        for contig_id, rest in embeddings.items():
            embeddings_contig_dict = embeddings[contig_id]

            for sequence_id, embedding in embeddings_contig_dict.items():
                # if proteins, don't make another dict
                if contig_id == "proteins":
                    embeddings_name = sequence_id
                else:
                    embeddings_name = f"{contig_id}:{sequence_id}"
                hf.create_dataset(embeddings_name, data=embedding)


def write_predictions(
    predictions: Dict[str, Dict[str, Tuple[List[str], Any, Any]]],
    out_path: Path,
    proteins_flag: bool,
    mask_threshold: float,
) -> None:
    """
    Write predictions to an output file.

    Args:
        predictions (Dict[str, Dict[str, Tuple[List[str], Any, Any]]]): Predictions dictionary containing contig IDs, sequence IDs, predictions, and additional information.
        out_path (Path): Path to the output file.
        proteins_flag (bool): Flag indicating whether the predictions are in proteins mode or not.
        mask_threshold (float): between 0 and 100 - below this ProstT5 confidence, 3Di predictions are masked


    Returns:
        None
    """
    ss_mapping = {
        0: "A",
        1: "C",
        2: "D",
        3: "E",
        4: "F",
        5: "G",
        6: "H",
        7: "I",
        8: "K",
        9: "L",
        10: "M",
        11: "N",
        12: "P",
        13: "Q",
        14: "R",
        15: "S",
        16: "T",
        17: "V",
        18: "W",
        19: "Y",
        20: "X",  # fully mask the low confidence 3Di residues with X not lower case (not working for Foldseek v10, but X does)
        # 20: "a",
        # 21: "c",
        # 22: "d",
        # 23: "e",
        # 24: "f",
        # 25: "g",
        # 26: "h",
        # 27: "i",
        # 28: "k",
        # 29: "l",
        # 30: "m",
        # 31: "n",
        # 32: "p",
        # 33: "q",
        # 34: "r",
        # 35: "s",
        # 36: "t",
        # 37: "v",
        # 38: "w",
        # 39: "y"
    }

    mask_prop_threshold = mask_threshold / 100

    with open(out_path, "w+") as out_f:
        for contig_id, rest in predictions.items():
            prediction_contig_dict = predictions[contig_id]

            # Filter out entries where the length of the value is 0
            # Issue #47
            prediction_contig_dict = {
                k: v for k, v in prediction_contig_dict.items() if len(v[0]) > 0
            }

            # masking - make the 3Di X=20
            for key, (pred, mean_prob, all_prob) in prediction_contig_dict.items():
                for i in range(len(pred)):
                    if all_prob[i] < mask_prop_threshold: # flat (L,)
                        pred[i] = 20

            if proteins_flag is True:
                # no contig_id
                out_f.write(
                    "".join(
                        [
                            ">{}\n{}\n".format(
                                f"{seq_id}",
                                "".join(
                                    list(map(lambda yhat: ss_mapping[int(yhat)], yhats))
                                ),
                            )
                            for seq_id, (yhats, _, _) in prediction_contig_dict.items()
                        ]
                    )
                )

            else:
                # writes each CDS to a 3di FASTA with header contig_id:CDS_id
                out_f.write(
                    "".join(
                        [
                            ">{}\n{}\n".format(
                                f"{contig_id}:{seq_id}",
                                "".join(
                                    list(map(lambda yhat: ss_mapping[int(yhat)], yhats))
                                ),
                            )
                            for seq_id, (yhats, _, _) in prediction_contig_dict.items()
                        ]
                    )
                )

    logger.info(f"Finished writing results to {out_path}")
    return None


def write_probs(
    predictions: Dict[str, Tuple[int, float, Union[int, np.ndarray]]],
    output_path_mean: Path,
    output_path_all: Path,
) -> None:
    """
    Write all ProstT5 encoder + CNN probabilities and mean probabilities to output files.

    Args:
        predictions (Dict[str, Tuple[int, float, Union[int, np.ndarray]]]):
        Predictions dictionary containing  sequence IDs, probabilities, and additional information.
        output_path_mean (str): Path to the output file for mean probabilities.
        output_path_all (str): Path to the output file for all probabilities.

    Returns:
        None
    """

    for contig_id, contig_predictions in predictions.items():
        # contig_predictions: Dict[seq_id, (labels, mean_prob, all_probs)]

        # ---- write mean probabilities ----
        with open(output_path_mean, "w+") as out_f:
            for seq_id, (N, mean_prob, N) in contig_predictions.items():

                out_f.write(f"{seq_id},{mean_prob}\n")

        # ---- write per-residue probabilities ----
        if output_path_all is not None:
            with open(output_path_all, "w+") as out_f:
                for seq_id, (N, N, all_probs) in contig_predictions.items():

                    # convert to percentage - no need to flatten as all_probs # flat (L,)
                    probs = (all_probs * 100).tolist()

                    rounded_probs = [round(p, 2) for p in probs]

                    out_f.write(
                        json.dumps(
                            {
                                "seq_id": seq_id,
                                "probability": rounded_probs,
                            }
                        )
                        + "\n"
                    )


def toCPU(tensor: torch.Tensor) -> np.ndarray:
    """
    Move a tensor to CPU and convert it to a NumPy array.

    Args:
        tensor (torch.Tensor): Input tensor.

    Returns:
        np.ndarray: NumPy array.
    """
    if len(tensor.shape) > 1:
        return tensor.detach().cpu().squeeze(dim=-1).numpy()
    else:
        return tensor.detach().cpu().numpy()


def load_predictor(checkpoint_path: Union[str, Path]) -> CNN:
    """
    Load a pre-trained CNN ProstT5 prediction head weights from a checkpoint file.

    Args:
        checkpoint_path (Union[str, Path]): Path to the checkpoint file.

    Returns:
        CNN: Loaded CNN model.
    """

    model = CNN()

    # checkpoint_path = Path(CNN_DIR) / "cnn_chkpnt" / "model.pt"

    state = torch.load(checkpoint_path, map_location=device)

    # regular ProstT5 CNN
    if checkpoint_path.suffix == ".pt":
        model.load_state_dict(state["state_dict"])
    # finetuned
    else:
        model.load_state_dict(state)

    model = model.eval()
    model = model.to(device)

    return model

def chunk_sequence(seq, max_len):
    """
    Yield (start, subseq) splitting seq into `max_len` nearly equal chunks.
    No overlap. Preserves order.
    """
    L = len(seq)
    n_chunks = math.ceil(L / max_len)
    base = L // n_chunks
    remainder = L % n_chunks

    start = 0
    for i in range(n_chunks):
        # distribute the extra 1 from remainder to the first `remainder` chunks
        size = base + (1 if i < remainder else 0)
        end = start + size
        yield start, seq[start:end]
        start = end


def get_embeddings(
    cds_dict: Dict[str, Dict[str, Tuple[str, ...]]],
    out_path: Path,
    prefix: str,
    model_dir: Path,
    model_name: str,
    checkpoint_path: Path,
    output_3di: Path,
    output_h5_per_residue: Path,
    output_h5_per_protein: Path,
    half_precision: bool,
    max_batch_residues: int = 100000,
    max_seq_len: int = 30000,
    max_batch: int = 10000,
    cpu: bool = False,
    output_probs: bool = True,
    proteins_flag: bool = False,
    save_per_residue_embeddings: bool = False,
    save_per_protein_embeddings: bool = False,
    threads: int = 1,
    mask_threshold: float = 0,
) -> bool:
    """
    Generate embeddings and predictions for protein sequences using ProstT5 encoder & CNN prediction head.

    Args:
        cds_dict (Dict[str, Dict[str, Tuple[str, ...]]]): nested dictionary containing contig IDs, CDS IDs and corresponding protein sequences.
        out_path (Path): Path to the output directory.
        prefix (str): Prefix for the output files.
        model_dir (Path): Directory containing the pre-trained model.
        model_name (str): Name of the pre-trained model.
        output_3di (Path): Path to the output 3Di file.
        output_h5_per_residue (Path): Path to the output h5 per residue embeddings file.
        output_h5_per_protein (Path): Path to the output h5 per proteins embeddings file.
        half_precision (bool): Whether to use half precision for the models.
        max_batch_residues (int, optional): Maximum number of residues allowed in a batch. 
        max_seq_len (int, optional): Maximum sequence length allowed. 
        max_batch (int, optional): Maximum batch size. 
        cpu (bool, optional): Whether to use CPU for processing. Defaults to False.
        output_probs (bool, optional): Whether to output probabilities. Defaults to True.
        proteins_flag (bool, optional): Whether the sequences are proteins. Defaults to False.
        save_embeddings (bool, optional): Whether to save embeddings to h5 file. Defaults to False. Will  save per residue embeddings
        per_protein_embeddings (bool, optional): Whether to save per protein mean embeddings to h5 file. Defaults to False.
        threads (int): number of cpu threads
        mask_threshold (float) : 0-100 - below this ProstT5 confidence threshold, these residues are masked


    Returns:
        bool: True if embeddings and predictions are generated successfully.
    """

    chunk_len = 1568 # hardcode chunk length of 1568 for inference

    predictions = {}

    if save_per_residue_embeddings:
        embeddings_per_residue = {}
    if save_per_protein_embeddings:
        embeddings_per_protein = {}

    prostt5_prefix = "<AA2fold>"

    model, vocab = get_model(model_dir, model_name, cpu, threads)
    if model_name != "gbouras13/modernprost-base":
        predictor = load_predictor(checkpoint_path)

    logger.info(f"Beginning {model_name} predictions")

    # modernprost always loaded in HP
    if model_name != "gbouras13/modernprost-base":
        if half_precision:
            model = model.half()
            predictor = predictor.half()
            logger.info("Using models in half-precision")
        else:
            logger.info("Using models in full-precision")


    # loop over each record in the cds dict
    fail_ids = []
    for record_id, seq_record_dict in cds_dict.items():
        # instantiate the nested dict
        predictions[record_id] = {}
        batch_predictions = {}
        # for modernprost chunking
        chunk_store = defaultdict(dict)
        chunk_embeddings_per_residue = defaultdict(dict)
        chunk_embeddings_per_protein = defaultdict(dict)

        # embeddings
        if save_per_residue_embeddings:
            batch_embeddings_per_residue = {}
        if save_per_protein_embeddings:
            batch_embeddings_per_protein = {}

        seq_dict = []
        for k, feat in seq_record_dict.items():
            v = feat.qualifiers.get("translation")
            if v and isinstance(v, str):
                seq = v.replace("U", "X").replace("Z", "X").replace("O", "X")
                seq_dict.append((k, seq, len(seq)))
            else:
                logger.warning(f"Protein header {k} is corrupt. It will be saved in fails.tsv")
                fail_ids.append(k)

        original_keys = list(seq_record_dict.keys())
        # --- sort once ---
        seq_dict.sort(key=lambda x: x[2], reverse=True)


        batch = list()
        n_res_batch = 0
        for seq_idx, (pdb_id, seq, slen) in enumerate(tqdm(seq_dict, desc=f"Predicting 3Di for {record_id}"), 1):


            # replace non-standard AAs
            seq = seq.replace("U", "X").replace("Z", "X").replace("O", "X")
            seq_len = len(seq)

            # only needed for old models 

            if model_name != "gbouras13/modernprost-base":

                seq = prostt5_prefix + " " + " ".join(list(seq))

                batch.append((pdb_id, seq, seq_len))
            
            
            # only do chunking for modernprost as dont want to deal with prefix char plus for ProstT5 and anyway we are swapping to modernprost
            else:

                if slen > chunk_len:
                    for chunk_idx, (start, subseq) in enumerate(chunk_sequence(seq, chunk_len)): # hardcoded to 1568
                        chunk_pid = f"{pdb_id}__chunk{chunk_idx}"
                        batch.append((chunk_pid, subseq, len(subseq)))
                        n_res_batch += len(subseq)
                else:
                    batch.append((pdb_id, seq, slen))


            # count residues in current batch and add the last sequence length to
            # avoid that batches with (n_res_batch > max_residues) get processed
            n_res_batch += slen

            if (
                len(batch) >= max_batch
                or n_res_batch >= max_batch_residues
                or seq_idx == len(seq_dict)
                or seq_len > max_seq_len
            ):
                pdb_ids, seqs, seq_lens = zip(*batch)
                batch.clear()
                n_res_batch = 0
                batch = list()

                # this is the same for modernprost and ProstT5
                token_encoding = vocab.batch_encode_plus(
                    seqs,
                    add_special_tokens=True,
                    padding="longest",
                    return_tensors="pt",
                ).to(device)

                try:

                    if model_name != "gbouras13/modernprost-base":


                        with torch.no_grad():
                            embedding_repr = model(
                                token_encoding.input_ids,
                                attention_mask=token_encoding.attention_mask,
                            )
                    else:
                        with torch.no_grad():
                            outputs = model(**token_encoding)


                except RuntimeError:
                    logger.warning(f" number of residues in batch {n_res_batch}")
                    logger.warning(f" seq length is {seq_len}")
                    logger.warning(f" ids are {pdb_ids}")
                    logger.warning(
                        "RuntimeError during embedding for {} (L={})".format(
                            pdb_id, seq_len
                        )
                    )
                    for id in pdb_ids:
                        fail_ids.append(id)
                    continue

                try:

                    if model_name != "gbouras13/modernprost-base":

                        # ProtT5 appends a special tokens at the end of each sequence
                        # Mask this also out during inference while taking into account the prostt5 prefix


                        for idx, s_len in enumerate(seq_lens):
                            token_encoding.attention_mask[idx, s_len + 1] = 0

                        # extract last hidden states (=embeddings)
                        residue_embedding = embedding_repr.last_hidden_state.detach()
                        # mask out padded elements in the attention output (can be non-zero) for further processing/prediction
                        residue_embedding = (
                            residue_embedding
                            * token_encoding.attention_mask.unsqueeze(dim=-1)
                        )
                        # slice off embedding of special token prepended before to each sequence
                        residue_embedding = residue_embedding[:, 1:]
                        prediction = predictor(residue_embedding)

                        # compute max probabilities per token/residue
                        probabilities = toCPU(
                            torch.max(F.softmax(prediction, dim=1), dim=1, keepdim=True)[0]
                        )

                        prediction = toCPU(
                            torch.max(prediction, dim=1, keepdim=True)[1]
                        ).astype(np.byte)


                        # batch-size x seq_len x embedding_dim
                        # extra token is added at the end of the seq
                        for batch_idx, identifier in enumerate(pdb_ids):
                            s_len = seq_lens[batch_idx]

                            # save embeddings
                            if save_per_residue_embeddings or save_per_protein_embeddings:
                                try:
                                    # account for prefix in offset
                                    emb = embedding_repr.last_hidden_state[
                                        batch_idx, 1 : s_len + 1
                                    ]

                                    if save_per_residue_embeddings:
                                        batch_embeddings_per_residue[identifier] = (
                                            emb.detach().cpu().numpy().squeeze()
                                        )

                                    if save_per_protein_embeddings:
                                        batch_embeddings_per_protein[identifier] = (
                                            emb.mean(dim=0).detach().cpu().numpy().squeeze()
                                        )

                                except:
                                    logger.warning(
                                        f"Saving embeddings failed for {identifier}"
                                    )

                            # slice off padding and special token appended to the end of the sequence
                            pred = prediction[batch_idx, :, 0:s_len].squeeze() # flat (L,) vector 

                            # always return the mean probs
                            mean_prob = round(
                                100 * np.mean(probabilities[batch_idx, :, 0:s_len]), 2
                            )

                            if output_probs:  # if you want the per-residue probs
                                all_prob = probabilities[batch_idx, :, 0:s_len].squeeze(0) #takes shape (1, L) and flattens to flat (L,) vector 
                                batch_predictions[identifier] = (
                                    pred,
                                    mean_prob,
                                    all_prob,
                                )
                            else:
                                batch_predictions[identifier] = (pred, mean_prob, None)

                            try:
                                len(batch_predictions[identifier][0])
                            except:
                                logger.warning(
                                    f"{identifier} {record_id} prediction has length 0"
                                )
                                fail_ids.append(identifier)
                                continue

                            if s_len != len(batch_predictions[identifier][0]):
                                logger.warning(
                                    f"Length mismatch for {identifier}: is:{len(batch_predictions[identifier][0])} vs should:{s_len}"
                                )

                    else:

                        # modernprost

                        logits = outputs.logits

                        probabilities = toCPU(
                            torch.max(F.softmax(logits, dim=-1), dim=-1, keepdim=True).values
                        )
                        for batch_idx, identifier in enumerate(pdb_ids): # way easier than ProstT5
                            s_len = seq_lens[batch_idx]


                            # save embeddings
                            if save_per_residue_embeddings or save_per_protein_embeddings:
                                try:
                                    # confusingly I called hidden_states = outputs.last_hidden_state https://huggingface.co/gbouras13/modernprost-base/blob/main/modeling_modernprost.py
                                    # so this is the last layer
                                    emb = outputs.hidden_states[
                                        batch_idx, 0:s_len
                                    ]

                                    if save_per_residue_embeddings:

                                        per_residue_embeddings = (emb.detach().cpu().numpy().squeeze())
                                        if "__chunk" in identifier:
                                            base_id, chunk_tag = identifier.split("__chunk")
                                            chunk_idx = int(chunk_tag)
                                            chunk_embeddings_per_residue[base_id][chunk_idx] = per_residue_embeddings
                                            
                                        else:

                                            batch_embeddings_per_residue[identifier] = per_residue_embeddings


                                    if save_per_protein_embeddings:
                                        per_protein_embeddings = (emb.mean(dim=0).detach().cpu().numpy().squeeze())
                                        if "__chunk" in identifier:
                                            base_id, chunk_tag = identifier.split("__chunk")
                                            chunk_idx = int(chunk_tag)
                                            chunk_embeddings_per_protein[base_id][chunk_idx] = per_protein_embeddings
                                            
                                        else:

                                            batch_embeddings_per_protein[identifier] = per_protein_embeddings


                                except:
                                    logger.warning(
                                        f"Saving embeddings failed for {identifier}"
                                    )



                            pred = logits[batch_idx, 0:s_len, :].squeeze()  # flat (L,) vector 

                            pred = toCPU(
                                torch.argmax(pred, dim=1, keepdim=True)
                            ).astype(np.byte)

                            all_prob = probabilities[batch_idx, 0:s_len]  # flat (L,) vector as probabilities is

                            if "__chunk" in identifier:
                                base_id, chunk_tag = identifier.split("__chunk")
                                chunk_idx = int(chunk_tag)

                                chunk_store[base_id][chunk_idx] = (
                                        pred,
                                        all_prob)

                                try:
                                    len(chunk_store[base_id][chunk_idx][0])
                                except:
                                    logger.warning(
                                        f"{identifier} {record_id} prediction has length 0"
                                    )
                                    fail_ids.append(base_id)
                                    continue

                                if s_len != len(chunk_store[base_id][chunk_idx][0]):
                                    logger.warning(
                                        f"Length mismatch for {identifier}: is:{len(chunk_store[base_id][chunk_idx][0])} vs should:{s_len}"
                                    )
                                
                            else:
                                mean_prob = round(100 * probabilities[batch_idx, 0:s_len].mean().item(), 2)
                           
                                if output_probs:  # if you want the per-residue probs
                                    batch_predictions[identifier] = (
                                        pred,
                                        mean_prob,
                                        all_prob,
                                    )
                                else:
                                    batch_predictions[identifier] = (pred, mean_prob, None)


                                # some catch statements regardless of model

                                try:
                                    len(batch_predictions[identifier][0])
                                except:
                                    logger.warning(
                                        f"{identifier} {record_id} prediction has length 0"
                                    )
                                    fail_ids.append(identifier)
                                    continue

                                if s_len != len(batch_predictions[identifier][0]):
                                    logger.warning(
                                        f"Length mismatch for {identifier}: is:{len(batch_predictions[identifier][0])} vs should:{s_len}"
                                    )


                    ######################################
                    # --- recombine chunked sequences ---

                    if model_name == "gbouras13/modernprost-base":
                        for pid, chunks in chunk_store.items():
                            preds = []
                            probs_all = []

                            for idx in sorted(chunks):
                                
                                pred, prob = chunks[idx]
                                preds.append(pred)
                                if prob is not None:
                                    probs_all.append(prob)

                            # always needed to calculate the mean confidence
                            pred_full = np.concatenate(preds)
                            probs_full = np.concatenate(probs_all)
                            mean_prob = round(100 * probs_full.mean(), 2)


                            if not output_probs: # only output full probs if true
                                probs_full = None

                            batch_predictions[pid] = (
                                pred_full,
                                mean_prob,
                                probs_full,
                            )
                      
                        if save_per_residue_embeddings:
                                        
                            for pid, chunks in chunk_embeddings_per_residue.items():
                                embs_all = []

                                for idx in sorted(chunks):
                                    emb = chunks[idx]
                                    embs_all.append(emb)

                                # concatenate the actual collected embeddings
                                embs_full = np.concatenate(embs_all, axis=0)  # make sure axis is correct

                                batch_embeddings_per_residue[pid] = embs_full

                        if save_per_protein_embeddings:
                                        
                            for pid, chunks in chunk_embeddings_per_protein.items():
                                embs_all = []

                                for idx in sorted(chunks):
                                    emb = chunks[idx]
                                    embs_all.append(emb)

                                # concatenate the actual collected embeddings
                                embs_full = np.concatenate(embs_all, axis=0)

                                batch_embeddings_per_protein[pid] = embs_full


                    # reorder to match the original FASTA
                    predictions[record_id] = {}

                    for k in original_keys:
                        if k in batch_predictions:
                            predictions[record_id][k] = batch_predictions[k]
                    
                    if save_per_residue_embeddings:
                        embeddings_per_residue[record_id] = {}
                        for k in original_keys:
                            if k in batch_predictions:
                                embeddings_per_residue[record_id][k] = batch_embeddings_per_residue[k]

                    if save_per_protein_embeddings:
                        embeddings_per_protein[record_id] = {}
                        for k in original_keys:
                            if k in batch_predictions:
                                embeddings_per_protein[record_id][k] = batch_embeddings_per_protein[k]


                except IndexError:
                    logger.warning(
                        "Index error during prediction for {} (L={})".format(
                            pdb_id, seq_len
                        )
                    )
                    for id in pdb_ids:
                        fail_ids.append(id)
                    continue

    # write list of fails if length > 0
    if len(fail_ids) > 0:
        fail_tsv: Path = Path(out_path) / "fails.tsv"

        # Convert the list to a list of lists
        data_as_list_of_lists = [[str(item)] for item in fail_ids]

        # Write the list to a TSV file
        with open(fail_tsv, "w", newline="") as file:
            tsv_writer = csv.writer(file, delimiter="\t")
            tsv_writer.writerows(data_as_list_of_lists)

    write_predictions(predictions, output_3di, proteins_flag, mask_threshold)

    if save_per_residue_embeddings:
        write_embeddings(embeddings_per_residue, output_h5_per_residue)

    if save_per_protein_embeddings:
        write_embeddings(embeddings_per_protein, output_h5_per_protein)

    # always write the mean embeddings
    if model_name == "gbouras13/modernprost-base":
        model_prefix = "modernprost"
    else:
        model_prefix = "prostT5"
    mean_probs_out_path: Path = (
        Path(out_path) / f"{prefix}_{model_prefix}_3di_mean_probabilities.csv"
    )

    # output per residue probs
    if output_probs:
        all_probs_out_path: Path = (
            Path(out_path) / f"{prefix}_{model_prefix}_3di_all_probabilities.json"
        )
    else:
        all_probs_out_path = None

    write_probs(predictions, mean_probs_out_path, all_probs_out_path, original_keys)

    return predictions

