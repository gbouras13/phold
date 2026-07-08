#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3Di prediction for phold — wraps pholdlib's shared inference engine.

Code adapted from @mheinzinger
https://github.com/mheinzinger/ProstT5/blob/main/scripts/predict_3Di_encoderOnly.py
"""

import csv
import json
from contextlib import ExitStack
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import h5py
import numpy as np

from loguru import logger

from phold.utils.util import atomic_write_path

# ── pholdlib shared components ────────────────────────────────────────────────
from pholdlib.prostt5.model import CNN, get_T5_model, load_predictor, toCPU  # noqa: F401
from pholdlib.prostt5.inference import (
    run_prostt5_inference,
    run_prostt5_inference_multi_gpu,
)
from pholdlib.prostt5.device import parse_gpus
from pholdlib.prostt5.output import (  # noqa: F401
    SS_MAPPING,
    mask_low_confidence_aa,
    write_fail_ids,
)

# ── phold-specific DB helpers ─────────────────────────────────────────────────
from phold.databases.db import check_prostT5_download, download_zenodo_prostT5


# ─────────────────────────────────────────────────────────────────────────────
# HDF5 embedding writers (phold-specific: nested contig→seq structure)
# ─────────────────────────────────────────────────────────────────────────────

def write_embeddings(
    embeddings: Dict[str, Dict[str, Any]],
    out_path: Path,
) -> None:
    """Write per-residue or per-protein embeddings to HDF5.

    Keys are stored as ``contig_id:seq_id`` (or just ``seq_id`` in proteins mode).

    Atomic: the HDF5 is written to a sibling temp file and only renamed onto
    ``out_path`` on success. Crashes mid-write (OOM, Ctrl-C, disk full) used to
    leave a truncated ``.h5`` at the final path; ``h5py`` doesn't validate the
    superblock when re-opening, so the next run would attempt to keep using
    the broken file. Now: temp gets unlinked on any ``BaseException`` and the
    original ``out_path`` (if any) is left byte-identical.
    """
    with atomic_write_path(out_path) as tmp:
        with h5py.File(str(tmp), "w") as hf:
            for contig_id, contig_dict in embeddings.items():
                for sequence_id, embedding in contig_dict.items():
                    key = sequence_id if contig_id == "proteins" else f"{contig_id}:{sequence_id}"
                    hf.create_dataset(key, data=embedding)


# ─────────────────────────────────────────────────────────────────────────────
# 3Di FASTA writer (phold-specific: nested contig structure + masking)
# ─────────────────────────────────────────────────────────────────────────────

def write_predictions(
    predictions: Dict[str, Dict[str, Tuple]],
    out_path: Path,
    proteins_flag: bool,
    mask_threshold: float,
) -> None:
    """Write 3Di predictions to a FASTA file.

    Args:
        predictions: Nested dict ``{contig_id: {seq_id: (pred, mean_prob, all_prob)}}``.
        out_path: Output FASTA path.
        proteins_flag: If True, headers are plain ``seq_id``; otherwise ``contig_id:seq_id``.
        mask_threshold: Residues whose max softmax probability (0–100) is below this
                        threshold are replaced with 'X'.
    """
    mask_prop = mask_threshold / 100

    with open(out_path, "w+") as out_f:
        for contig_id, contig_dict in predictions.items():
            # drop zero-length predictions (issue #47)
            contig_dict = {k: v for k, v in contig_dict.items() if len(v[0]) > 0}

            # apply masking in-place: pred is np.byte, all_prob shape (1, L)
            for key, (pred, mean_prob, all_prob) in contig_dict.items():
                pred[all_prob[0] < mask_prop] = 20  # 'X'

            header_fmt = "{}" if proteins_flag else "{}:{}"

            out_f.write(
                "".join(
                    ">{}\n{}\n".format(
                        header_fmt.format(seq_id) if proteins_flag
                        else header_fmt.format(contig_id, seq_id),
                        "".join(SS_MAPPING[int(y)] for y in yhats),
                    )
                    for seq_id, (yhats, _, _) in contig_dict.items()
                )
            )

    logger.info(f"Finished writing 3Di FASTA to {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Probability writers — bioconda-compatible signature (takes cds_dict)
# ─────────────────────────────────────────────────────────────────────────────

def write_probs(
    predictions: Dict[str, Dict[str, Tuple]],
    output_path_mean: Path,
    output_path_all: Optional[Path],
    cds_dict: Dict[str, Dict],
) -> None:
    """Write mean and per-residue probabilities to disk.

    The file is opened once and all contigs are written sequentially so that
    multi-contig GenBank inputs produce a single combined CSV / JSONL.

    Args:
        predictions: Nested ``{contig_id: {seq_id: (labels, mean_prob, all_probs)}}``.
        output_path_mean: CSV path — one ``seq_id,mean_prob`` line per sequence.
        output_path_all:  JSONL path for per-residue probs. Pass ``None`` to skip.
        cds_dict: The original CDS dict; used to retrieve the per-contig key order.
    """
    # Atomicity + leak safety in one pass:
    #
    # * Previously ``out_all`` was opened OUTSIDE the ``with``; if any
    #   ``out_mean.write`` or ``out_all.write`` raised mid-loop, the second
    #   handle leaked and the manual ``close()`` at the end never ran.
    # * Both files were also written directly to their final paths. A
    #   crash mid-loop left a partially-written CSV / JSONL on disk, and
    #   ``--restart`` happily picked it up as if it were complete.
    #
    # ``ExitStack`` registers both atomic-write contexts and their inner
    # file handles for LIFO cleanup, so on any failure (incl. KeyboardInterrupt
    # — atomic_write_path catches BaseException) the file handles close, the
    # sibling temp files are unlinked, and the final paths are left exactly
    # as they were before the call. On success, the file handles close
    # first (flushing buffers) and then atomic_write_path's ``__exit__``
    # ``os.replace``s the temp onto the final path — atomic for any
    # downstream observer.
    with ExitStack() as stack:
        tmp_mean = stack.enter_context(atomic_write_path(output_path_mean))
        out_mean = stack.enter_context(open(tmp_mean, "w"))

        if output_path_all is not None:
            tmp_all = stack.enter_context(atomic_write_path(output_path_all))
            out_all = stack.enter_context(open(tmp_all, "w"))
        else:
            out_all = None

        for contig_id, contig_predictions in predictions.items():
            original_keys = list(cds_dict[contig_id].keys())

            for seq_id in original_keys:
                if seq_id not in contig_predictions:
                    logger.warning(f"Missing ProstT5 mean confidence for {seq_id}")
                    continue

                _, mean_prob, all_probs = contig_predictions[seq_id]
                # NOTE: mean_prob is np.float32 (see inference.py TODO). The f-string
                # deliberately formats it as the full float32 repr (e.g. 82.69000244140625)
                # to match bioconda 1.2.5 byte-for-byte. Once the inference.py TODO is
                # resolved and mean_prob becomes a clean Python float, this line will
                # automatically produce clean 2dp output (e.g. 82.69) — no change needed here.
                out_mean.write(f"{seq_id},{mean_prob}\n")

                if out_all is not None:
                    # TODO: .tolist() converts float32 → Python float (float64), then
                    # round(p, 2) gives clean 2dp Python floats. Verify this matches
                    # bioconda 1.2.5's all_probabilities.json output once a run completes.
                    probs = (all_probs * 100).flatten().tolist()
                    rounded_probs = [round(p, 2) for p in probs]
                    out_all.write(
                        json.dumps({"seq_id": seq_id, "probability": rounded_probs}) + "\n"
                    )


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def get_embeddings(
    cds_dict: Dict[str, Dict[str, Any]],
    out_path: Path,
    prefix: str,
    model_dir: Path,
    model_name: str,
    checkpoint_path: Path,
    output_3di: Path,
    output_h5_per_residue: Path,
    output_h5_per_protein: Path,
    half_precision: bool,
    max_residues: int = 100000,
    max_seq_len: int = 30000,
    max_batch: int = 10000,
    cpu: bool = False,
    output_probs: bool = True,
    proteins_flag: bool = False,
    save_per_residue_embeddings: bool = False,
    save_per_protein_embeddings: bool = False,
    threads: int = 1,
    mask_threshold: float = 0,
    gpus: Optional[str] = None,
) -> Dict:
    """Run ProstT5 + CNN 3Di prediction for all sequences in *cds_dict*.

    Args:
        cds_dict: Nested ``{contig_id: {seq_id: BioPython_feature}}``.
                  Feature must expose ``qualifiers["translation"]``.
        out_path: Directory for output files.
        prefix: Filename prefix for CSV / JSONL outputs.
        model_dir: Directory where ProstT5 is cached.
        model_name: HuggingFace model identifier.
        checkpoint_path: Path to the CNN ``.pt`` checkpoint.
        output_3di: Output FASTA path for 3Di sequences.
        output_h5_per_residue: HDF5 path for per-residue embeddings.
        output_h5_per_protein: HDF5 path for per-protein embeddings.
        half_precision: If True, cast model + predictor to fp16 after loading.
        max_residues: Max total residues per inference batch.
        max_seq_len: Sequences longer than this flush a batch immediately.
        max_batch: Max sequences per batch.
        cpu: Force CPU inference.
        output_probs: Whether to write per-residue probability JSONL.
        proteins_flag: True when input is a flat proteins FASTA (no contig nesting).
        save_per_residue_embeddings: Save per-residue HDF5.
        save_per_protein_embeddings: Save per-protein HDF5.
        threads: Number of CPU threads for torch.
        mask_threshold: Residues with max softmax prob < threshold/100 → 'X'.
        gpus: Comma-separated CUDA indices (e.g. "0,2"). None = auto-detect
              all visible CUDA GPUs (or fall through to MPS/XPU/CPU).
              Overridden by ``cpu=True``.

    Returns:
        predictions: Nested ``{contig_id: {seq_id: (pred, mean_prob, all_prob)}}``.
    """
    # ── resolve devices once, here ──────────────────────────────────────────
    devices = parse_gpus(cpu, gpus)
    logger.info(f"Beginning ProstT5 predictions on device(s): {devices}")
    if half_precision and devices == ["cpu"]:
        logger.info("CPU device — forcing full-precision (half-precision disabled).")
        half_precision = False
    if half_precision:
        logger.info("Using models in half-precision")
    else:
        logger.info("Using models in full-precision")

    # ── flatten contigs into one sorted seq_dict, keep a reverse map ────────
    # Multi-GPU sharding works best when the sort is global so each shard
    # contains a mix of long+short sequences. We re-nest after inference.
    flat_seq_dict: List[Tuple[str, str, int]] = []
    seq_to_contig: Dict[str, str] = {}
    all_fail_ids: List[str] = []

    for record_id, seq_record_dict in cds_dict.items():
        for k, feat in seq_record_dict.items():
            v = feat.qualifiers.get("translation")
            if v and isinstance(v, str):
                flat_seq_dict.append((k, v, len(v)))
                seq_to_contig[k] = record_id
            else:
                logger.warning(
                    f"Protein header {k} is corrupt. It will be saved in fails.tsv"
                )
                all_fail_ids.append(k)

    # sort once globally — every shard then begins length-descending
    flat_seq_dict.sort(key=lambda x: x[2], reverse=True)

    # ── inference — single-device or multi-GPU, both go through dispatcher ──
    flat_preds, flat_emb_res, flat_emb_prot, fail_ids = run_prostt5_inference_multi_gpu(
        flat_seq_dict,
        devices=devices,
        model_dir=model_dir,
        model_name=model_name,
        checkpoint_path=checkpoint_path,
        half_precision=half_precision,
        threads=threads,
        check_fn=check_prostT5_download,
        zenodo_fn=download_zenodo_prostT5,
        max_residues=max_residues,
        max_seq_len=max_seq_len,
        max_batch=max_batch,
        output_probs=output_probs,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        desc="Predicting 3Di",
    )
    all_fail_ids.extend(fail_ids)

    # ── re-nest by contig and restore original per-contig key order ─────────
    predictions: Dict = {rid: {} for rid in cds_dict}
    embeddings_per_residue: Dict = {rid: {} for rid in cds_dict} if save_per_residue_embeddings else {}
    embeddings_per_protein: Dict = {rid: {} for rid in cds_dict} if save_per_protein_embeddings else {}

    for record_id, seq_record_dict in cds_dict.items():
        original_keys = list(seq_record_dict.keys())
        predictions[record_id] = {
            k: flat_preds[k] for k in original_keys if k in flat_preds
        }
        if save_per_residue_embeddings:
            embeddings_per_residue[record_id] = {
                k: flat_emb_res[k] for k in original_keys if k in flat_emb_res
            }
        if save_per_protein_embeddings:
            embeddings_per_protein[record_id] = {
                k: flat_emb_prot[k] for k in original_keys if k in flat_emb_prot
            }

    # ── write outputs ────────────────────────────────────────────────────────
    if all_fail_ids:
        fail_tsv = Path(out_path) / "fails.tsv"
        write_fail_ids(all_fail_ids, fail_tsv)

    write_predictions(predictions, output_3di, proteins_flag, mask_threshold)

    if save_per_residue_embeddings:
        write_embeddings(embeddings_per_residue, output_h5_per_residue)
    if save_per_protein_embeddings:
        write_embeddings(embeddings_per_protein, output_h5_per_protein)

    mean_probs_path = Path(out_path) / f"{prefix}_prostT5_3di_mean_probabilities.csv"
    all_probs_path = (
        Path(out_path) / f"{prefix}_prostT5_3di_all_probabilities.json"
        if output_probs else None
    )

    write_probs(predictions, mean_probs_path, all_probs_path, cds_dict)

    return predictions
