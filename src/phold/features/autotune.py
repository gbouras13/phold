
# autobatch

from phold.features.predict_3Di import  get_T5_model
from tqdm import tqdm
import random
import torch
import torch.nn.functional as F
from loguru import logger
import time
import math
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from phold.io.handle_genbank import open_protein_fasta_file

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

    # get device
    device = None

    if cpu is True:
        device = torch.device("cpu")
    else:
        # check for NVIDIA/cuda
        if torch.cuda.is_available():
            device = torch.device("cuda:0")
        # check for apple silicon/metal
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
        else:
            device = torch.device("cpu")


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
            for i in tqdm(range(0, len(probe_seqs), bs), desc="Processing"):
                batch_seqs = probe_seqs[i : i + bs]

                print(batch_seqs)
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

            logger.info(f"Time elapsed {round(total_time,5)}")
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

    logger.info(f"##########################")
    logger.info(f"Best batch size: {best_bs}")
    # logger.info(f"best max residues: {best_residues}")

    return best_bs, best_residues



def run_autotune(    
    input_path,
    model_dir,
    model_name,
    cpu,
    threads,
    step, 
    min_batch,
    max_batch, 
    sample_seqs):

    # Dictionary to store the records
    cds_dict = {}


    with open_protein_fasta_file(input_path) as handle:  # handles gzip too
        records = list(SeqIO.parse(handle, "fasta"))
        if not records:
            logger.warning(f"No proteins were found in your input file {input_path}.")
            logger.error(
                f"Your input file {input_path} is likely not a amino acid FASTA file. Please check this."
            )
        for record in records:
            prot_id = record.id
            feature_location = FeatureLocation(0, len(record.seq))
            # Seq needs to be saved as the first element in list hence the closed brackets [str(record.seq)]
            seq_feature = SeqFeature(
                feature_location,
                type="CDS",
                qualifiers={
                    "ID": record.id,
                    "description": record.description,
                    "translation": str(record.seq),
                },
            )

            cds_dict[prot_id] = seq_feature

    if not cds_dict:
        logger.error(f"Error: no AA protein sequences found in {input_path} file")


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
        start_bs=min_batch,
        max_bs=max_batch,
        step=step # step size
    )

    logger.info(f"Optimal batch size is {batch_size} (residues per batch {max_residues})")

    return batch_size
