
# autobatch

from predict_3Di import device, get_T5_model
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import random
import h5py
import numpy as np
import torch
import torch.nn.functional as F
from loguru import logger
from torch import nn
from tqdm import tqdm
from transformers import T5EncoderModel, T5Tokenizer
import time
from phold.databases.db import check_prostT5_download, download_zenodo_prostT5
from phold.utils.constants import CNN_DIR
import math


"""

all_phold_structures.fasta is the FASTA of all the Phold DB 1.36M sequences

seqkit sample -n 5000 -j 4 all_phold_structures.fasta > all_phold_structures_5000.fasta

seqkit stats all_phold_structures_5000.fasta 
file                             format  type     num_seqs    sum_len  min_len  avg_len  max_len
all_phold_structures_5000.fasta  FASTA   Protein     5,000  1,091,959       21    218.4    2,793

# seems reasonable

"""

def sample_probe_sequences(seqs, n=5000, seed=0):
    """
    samples sequences 
    
    """

    rng = random.Random(seed)

    if n >= len(seqs):
        sampled = list(seqs)
    else:
        sampled = rng.sample(seqs, n)

    # sort by sequence length
    sampled.sort(key=len, reverse=True)

    return sampled

def autotune_batching_real_data(
    model_dir,
    model_name,
    cpu,
    threads,
    probe_seqs,
    start_bs=1,
    max_bs=150,
    step=10 # step size
):
    
    model, tokenizer = get_T5_model(model_dir, model_name, cpu, threads)

    model.eval()
    model.half()

    bs = start_bs
    results = []

    while bs <= max_bs:
        try:
            
            # seqs = probe_seqs
            n_tokens = sum(len(s) for s in probe_seqs)

            logger.info(f"Running with batch size {bs}")

            model.eval()

            total_tokens = 0
            total_time = 0.0
            batches = 0

            # iterate over real sequences in batches
            for i in range(0, len(probe_seqs), bs):
                batch_seqs = probe_seqs[i : i + bs]
                n_tokens = sum(len(s) for s in batch_seqs)
                total_tokens += n_tokens

                inputs = tokenizer(
                    batch_seqs,
                    padding=True,
                    return_tensors="pt",
                )
                inputs.pop("token_type_ids", None)
                inputs = {k: v.to(device) for k, v in inputs.items()}

                # timing
                torch.cuda.synchronize()
                t0 = time.perf_counter()
                with torch.no_grad():
                    _ = model(**inputs)
                torch.cuda.synchronize()

                total_time += time.perf_counter() - t0
                
                batches += 1

            time_per_token = total_time / total_tokens


            token_per_batch = math.floor(total_tokens / batches)

        
            results.append({
                "bs": bs,
                "tokens_per_batch": token_per_batch,
                "time": total_time,
                "time_per_token": time_per_token,
            })

            logger.info(f"Time elapsed {total_time}")
            logger.info(f"Tokens per batch {token_per_batch}")

            bs += step

        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            break

    
    if not results:
        raise RuntimeError("No batch size fits on this GPU")

    best_entry = min(results, key=lambda x: x["time_per_token"])

    best_bs = best_entry["bs"]
    best_residues = best_entry["tokens_per_batch"]
    # best_tpt = best_bs["time_per_token"]

    logger.info(f"best batch size: {best_bs}")
    # logger.info(f"best max residues: {best_residues}")

    return best_bs, best_residues



def run_autotune(    model_dir,
    model_name,
    cpu,
    threads,
    cds_dict, step, max_batch, sample_seqs):


    seqs = []
    for feat in cds_dict.values():
        v = feat.qualifiers.get("translation")
        if v and isinstance(v, str):
            seqs.append(v)

    logger.info("Beginning batch size tuning")
    logger.info(f"Using minimum batch size of 1 and maximum batch size of {max_batch}")

    # define the sampling

    probe_seqs = sample_probe_sequences(seqs, n=sample_seqs)

    batch_size, max_residues = autotune_batching_real_data(
        model_dir,
        model_name,
        cpu,
        threads,
        probe_seqs,
        start_bs=1,
        max_bs=max_batch,
        step=step # step size
    )

    logger.info(f"Optimal batch size is {batch_size} (residues per batch {max_residues})")

    return batch_size
