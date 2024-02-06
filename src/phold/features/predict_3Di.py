#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Adapted from @mheinzinger 

https://github.com/mheinzinger/ProstT5/blob/main/scripts/predict_3Di_encoderOnly.py

"""

import csv
import json
import os
import shutil
import time
from pathlib import Path
from urllib import request

import numpy as np
import torch
import torch.nn.functional as F
from loguru import logger
from torch import nn
from transformers import T5EncoderModel, T5Tokenizer

from phold.utils.constants import MODEL_DB


# Convolutional neural network (two convolutional layers)
class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()

        self.classifier = nn.Sequential(
            nn.Conv2d(1024, 32, kernel_size=(7, 1), padding=(3, 0)),  # 7x32
            nn.ReLU(),
            nn.Dropout(0.0),
            nn.Conv2d(32, 20, kernel_size=(7, 1), padding=(3, 0)),
        )

    def forward(self, x):
        """
        L = protein length
        B = batch-size
        F = number of features (1024 for embeddings)
        N = number of classes (20 for 3Di)
        """
        x = x.permute(0, 2, 1).unsqueeze(
            dim=-1
        )  # IN: X = (B x L x F); OUT: (B x F x L, 1)
        Yhat = self.classifier(x)  # OUT: Yhat_consurf = (B x N x L x 1)
        Yhat = Yhat.squeeze(dim=-1)  # IN: (B x N x L x 1); OUT: ( B x L x N )
        return Yhat


def get_T5_model(model_dir, model_name, cpu, finetuned_model_path, finetune_flag):
    # sets the device

    # Torch load will map back to device from state, which often is GPU:0.
    # to overcome, need to explicitly map to active device

    global device

    if torch.cuda.is_available():
        if cpu is True:
            device = torch.device("cpu")
            dev_name = "cpu"
        else:
            device = torch.device("cuda:0")
            dev_name = "cuda:0"
    else:
        device = torch.device("cpu")
        dev_name = "cpu"

    # logger device only if the function is called
    logger.info("Using device: {}".format(dev_name))

    # make dir
    Path(model_dir).mkdir(parents=True, exist_ok=True)
    # set as cache dir
    # os.environ['TRANSFORMERS_CACHE'] = f"{model_dir}/"
    # load
    logger.info(f"Loading T5 from: {model_dir}/{model_name}")
    logger.info(f"If {model_dir}/{model_name} is not found, it will be downloaded.")
    model = T5EncoderModel.from_pretrained(model_name, cache_dir=f"{model_dir}/").to(
        device
    )

    # finetuned weights
    if finetune_flag is True:
        # Load the non-frozen parameters from the saved file
        non_frozen_params = torch.load(finetuned_model_path, map_location=device)

        # Assign the non-frozen parameters to the corresponding parameters of the model
        for param_name, param in model.named_parameters():
            if param_name in non_frozen_params:
                param.data = non_frozen_params[param_name].data

    model = model.eval()
    vocab = T5Tokenizer.from_pretrained(
        model_name, cache_dir=f"{model_dir}/", do_lower_case=False
    )

    return model, vocab


def write_predictions(predictions, out_path):
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
    }

    with open(out_path, "w+") as out_f:
        for contig_id, rest in predictions.items():
            prediction_contig_dict = predictions[contig_id]

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


def write_probs(predictions, out_path_mean, output_path_all):
    with open(out_path_mean, "w+") as out_f:
        for contig_id, rest in predictions.items():
            prediction_contig_dict = predictions[contig_id]

            out_f.write(
                "\n".join(
                    [
                        "{},{}".format(seq_id, mean_prob)
                        for seq_id, (N, mean_prob, N) in prediction_contig_dict.items()
                    ]
                )
            )

    with open(output_path_all, "w+") as out_f:
        for contig_id, rest in predictions.items():
            prediction_contig_dict = predictions[contig_id]

            for seq_id, (N, N, all_probs) in prediction_contig_dict.items():
                # * 100
                all_probs = all_probs * 100
                # Convert NumPy array to list
                all_probs_list = (
                    all_probs.flatten().tolist()
                    if isinstance(all_probs, np.ndarray)
                    else all_probs
                )

                # round to 2 dp
                rounded_list = [round(num, 2) for num in all_probs_list]

                # Create a dictionary for the specific items
                specific_data = {"seq_id": seq_id, "probability": rounded_list}

                # Convert the dictionary to a JSON string
                json_data = json.dumps(specific_data)

                # Write the JSON string to the file
                out_f.write(json_data + "\n")  # Add a newline after each JSON object

    return None


def toCPU(tensor):
    if len(tensor.shape) > 1:
        return tensor.detach().cpu().squeeze(dim=-1).numpy()
    else:
        return tensor.detach().cpu().numpy()


def download_file(url, local_path):
    if not local_path.parent.is_dir():
        local_path.parent.mkdir()

    print("Downloading: {}".format(url))
    req = request.Request(
        url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 6.1; Win64; x64)"}
    )

    with request.urlopen(req) as response, open(local_path, "wb") as outfile:
        shutil.copyfileobj(response, outfile)
    return None


def load_predictor(
    weights_link="https://rostlab.org/~deepppi/prostt5/cnn_chkpnt/model.pt",
):
    model = CNN()
    checkpoint_p = Path(MODEL_DB) / "cnn_chkpnt" / "model.pt"

    # link seems broken
    # # if no pre-trained model is available, yet --> download it
    # if not checkpoint_p.exists():
    #     download_file(weights_link, checkpoint_p)

    state = torch.load(checkpoint_p, map_location=device)

    model.load_state_dict(state["state_dict"])

    model = model.eval()
    model = model.to(device)

    return model


def get_embeddings(
    cds_dict: dict,
    out_path,
    model_dir: Path,
    model_name: str,
    half_precision: bool,
    max_residues: int = 3000,
    max_seq_len: int = 1000,
    max_batch: int = 100,
    proteins: bool = False,
    cpu: bool = False,
    output_probs: bool = True,
    finetune_flag: bool = True,
    finetuned_model_path: str = None,
) -> bool:
    predictions = {}

    prefix = "<AA2fold>"

    if finetuned_model_path is None:
        finetuned_model_path = Path(MODEL_DB) / "Phrostt5_finetuned.pth"

    model, vocab = get_T5_model(
        model_dir, model_name, cpu, finetuned_model_path, finetune_flag
    )
    predictor = load_predictor(model_dir)

    logger.info("Beginning ProstT5 predictions.")

    if half_precision:
        model = model.half()
        predictor = predictor.half()
        logger.info("Using models in half-precision.")
    else:
        logger.info("Using models in full-precision.")

    # loop over each record in the cds dict
    fail_ids = []
    for record_id, cds_records in cds_dict.items():
        # instantiate the nested dict
        predictions[record_id] = {}

        seq_record_dict = cds_dict[record_id]
        seq_dict = {}

        # gets the seq_dict with key for id and the translation
        for key, seq_feature in seq_record_dict.items():
            # get the protein seq for normal
            if proteins is False:
                seq_dict[key] = seq_feature.qualifiers["translation"][0]
            else:  # proteins mode - it will be already in a dictionary
                seq_dict = seq_record_dict

        # sort sequences by length to trigger OOM at the beginning

        seq_dict = dict(
            sorted(seq_dict.items(), key=lambda kv: len(kv[1][0]), reverse=True)
        )

        # print("Average sequence length: {}".format(avg_length))
        # print("Number of sequences >{}: {}".format(max_seq_len, n_long))

        start = time.time()
        batch = list()
        for seq_idx, (pdb_id, seq) in enumerate(seq_dict.items(), 1):
            # print(pdb_id)
            # print(seq)

            # replace non-standard AAs
            seq = seq.replace("U", "X").replace("Z", "X").replace("O", "X")
            seq_len = len(seq)
            seq = prefix + " " + " ".join(list(seq))
            batch.append((pdb_id, seq, seq_len))

            # count residues in current batch and add the last sequence length to
            # avoid that batches with (n_res_batch > max_residues) get processed
            n_res_batch = sum([s_len for _, _, s_len in batch]) + seq_len
            if (
                len(batch) >= max_batch
                or n_res_batch >= max_residues
                or seq_idx == len(seq_dict)
                or seq_len > max_seq_len
            ):
                pdb_ids, seqs, seq_lens = zip(*batch)
                batch = list()

                token_encoding = vocab.batch_encode_plus(
                    seqs,
                    add_special_tokens=True,
                    padding="longest",
                    return_tensors="pt",
                ).to(device)
                try:
                    with torch.no_grad():
                        embedding_repr = model(
                            token_encoding.input_ids,
                            attention_mask=token_encoding.attention_mask,
                        )
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

                # ProtT5 appends a special tokens at the end of each sequence
                # Mask this also out during inference while taking into account the prefix
                try:
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

                    if output_probs:
                        # compute max probabilities per token/residue if requested
                        probabilities = toCPU(
                            torch.max(
                                F.softmax(prediction, dim=1), dim=1, keepdim=True
                            )[0]
                        )

                    prediction = toCPU(
                        torch.max(prediction, dim=1, keepdim=True)[1]
                    ).astype(np.byte)

                    # batch-size x seq_len x embedding_dim
                    # extra token is added at the end of the seq
                    for batch_idx, identifier in enumerate(pdb_ids):
                        s_len = seq_lens[batch_idx]

                        # # slice off padding and special token appended to the end of the sequence
                        # predictions[record_id][identifier] = prediction[
                        #     batch_idx, :, 0:s_len
                        # ].squeeze()

                        # slice off padding and special token appended to the end of the sequence
                        pred = prediction[batch_idx, :, 0:s_len].squeeze()

                        if output_probs:  # average over per-residue max.-probabilities
                            mean_prob = round(
                                100 * np.mean(probabilities[batch_idx, :, 0:s_len]), 2
                            )
                            all_prob = probabilities[batch_idx, :, 0:s_len]
                            predictions[record_id][identifier] = (
                                pred,
                                mean_prob,
                                all_prob,
                            )
                        else:
                            predictions[record_id][identifier] = (pred, None, None)

                        try:
                            len(predictions[record_id][identifier][0])
                        except:
                            logger.warning(
                                f"{identifier} {record_id} prediction has length 0"
                            )
                            fail_ids.append(identifier)
                            continue

                        if s_len != len(predictions[record_id][identifier][0]):
                            logger.warning(
                                f"Length mismatch for {identifier}: is:{len(predictions[record_id][identifier][0])} vs should:{s_len}"
                            )

                except IndexError:
                    logger.warning(
                        "Index error during prediction for {} (L={})".format(
                            pdb_id, seq_len
                        )
                    )
                    for id in pdb_ids:
                        fail_ids.append(id)
                    continue

    output_3di: Path = Path(out_path) / "output3di.fasta"

    # write list of fails if length > 0
    if len(fail_ids) > 0:
        fail_tsv: Path = Path(out_path) / "fails.tsv"

        # Convert the list to a list of lists
        data_as_list_of_lists = [[str(item)] for item in fail_ids]

        # Write the list to a TSV file
        with open(fail_tsv, "w", newline="") as file:
            tsv_writer = csv.writer(file, delimiter="\t")
            tsv_writer.writerows(data_as_list_of_lists)

    write_predictions(predictions, output_3di)

    mean_probs_out_path: Path = Path(out_path) / "output3di_mean_probabilities.csv"
    all_probs_out_path: Path = Path(out_path) / "output3di_all_probabilities.json"

    if output_probs:
        write_probs(predictions, mean_probs_out_path, all_probs_out_path)

    return True
