#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ModernProst 3Di + 12-state prediction for phold — wraps pholdlib's engine.

Where :mod:`phold.features.predict_3Di` runs ProstT5 plus a CNN head and emits
3Di only, this module runs a ModernProst checkpoint that emits Foldseek 3Di
*and* a 12-state secondary-structure alphabet from a single forward pass.

Two tasks are supported, following the checkpoint's training:

``classification`` (``modernprost-base`` / ``modernprost-50M``)
    argmax both heads and write ``{prefix}_3di.fasta`` + ``{prefix}_12st.fasta``.
    ``phold compare`` packs those into one combined Foldseek ``_ss`` database.

``pssm`` (``modernprost-pssm`` / ``modernprost-50M-pssm``)
    keep the full per-residue distributions and write them as profile text.
    ``phold compare`` searches Foldseek profile databases built from them.
    The argmax FASTAs are still written, so the outputs stay inspectable and
    the per-CDS confidence column works the same way in both tasks.

Masking note: ``--mask_threshold`` cannot be applied to the 3Di string here.
The combined encoding packs ``3di_index * 12 + ss12_index`` into one byte and
has no 21st "masked" 3Di state, so a masked ``X`` is unrepresentable. Masking is
therefore applied to the amino-acid FASTA only (as it already is on the ProstT5
path), which is what actually feeds Foldseek's amino-acid channel.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from loguru import logger

# Only light imports at module level. `pholdlib.modernprost.inference` pulls in
# torch and `.model` pulls in transformers — seconds locally, minutes on a
# shared cluster filesystem. `phold.subcommands.compare` imports this module
# just for `mean_probs_filename`, and never runs inference, so the heavy
# imports live inside `get_modernprost_predictions` instead.
from pholdlib.databases.modernprost import (
    CLASSIFICATION_TASK,
    PSSM_TASK,
    resolve_modernprost_model,
)

from phold.utils.util import atomic_write_path


# Filename component distinguishing ModernProst outputs from ProstT5's, so a
# directory can hold both without one silently overwriting the other.
MODEL_PREFIX = "modernprost"


def mean_probs_filename(prefix: str) -> str:
    """Name of the per-CDS mean-confidence CSV written by this module."""
    return f"{prefix}_{MODEL_PREFIX}_3di_mean_probabilities.csv"


def get_modernprost_predictions(
    cds_dict: Dict[str, Dict[str, Any]],
    out_path: Path,
    prefix: str,
    model_dir: Path,
    model_name: str,
    output_3di: Path,
    output_12st: Path,
    output_h5_per_residue: Path,
    output_h5_per_protein: Path,
    half_precision: bool,
    task: str = CLASSIFICATION_TASK,
    max_residues: int = 50000,
    max_seq_len: int = 30000,
    max_batch: int = 10000,
    chunk_len: Optional[int] = None,
    cpu: bool = False,
    output_probs: bool = True,
    proteins_flag: bool = False,
    save_per_residue_embeddings: bool = False,
    save_per_protein_embeddings: bool = False,
    threads: int = 1,
    gpus: Optional[str] = None,
) -> Tuple[Dict, Dict, Dict]:
    """Run ModernProst over every CDS in *cds_dict* and write phold's outputs.

    Args:
        cds_dict: Nested ``{contig_id: {seq_id: BioPython_feature}}``. Each
            feature must expose ``qualifiers["translation"]`` as a string.
        out_path: Directory for output files.
        prefix: Filename prefix for CSV / JSON outputs.
        model_dir: Directory where the model is cached.
        model_name: ModernProst registry name or full HuggingFace id.
        output_3di / output_12st: Output FASTA paths.
        output_h5_per_residue / output_h5_per_protein: HDF5 embedding paths.
        half_precision: Cast the model to fp16 (ignored on CPU).
        task: ``"classification"`` or ``"pssm"``.
        max_residues: Max total residues per inference batch.
        max_seq_len: Chunks longer than this flush a batch immediately.
        max_batch: Max sequences per batch.
        chunk_len: Sequences longer than this are split and reassembled.
            None uses pholdlib's default. It cannot be the default argument
            value here: reading it would import torch at module-import time,
            which is exactly what this module defers.
        cpu: Force CPU inference.
        output_probs: Whether to write the per-residue probability JSON.
        proteins_flag: True when input is a flat proteins FASTA (no contigs).
        save_per_residue_embeddings / save_per_protein_embeddings: HDF5 dumps.
        threads: Number of CPU threads for torch.
        gpus: Comma-separated CUDA indices. None = auto-detect.

    Returns:
        ``(predictions_3di, predictions_12st, profiles)`` — the first two nested
        ``{contig_id: {seq_id: (pred, mean_prob, all_prob)}}``, and ``profiles``
        a nested ``{contig_id: {seq_id: {"3di": ndarray, "12st": ndarray}}}``
        that is empty unless ``task="pssm"``.
    """
    # Deferred: torch / transformers (see the module docstring's note).
    from pholdlib.modernprost.inference import (
        DEFAULT_CHUNK_LEN, ModernProstResult, predictions_from_profiles,
        run_modernprost_inference_multi_gpu)
    from pholdlib.modernprost.output import (SS12_ALPHABET, THREEDI_ALPHABET,
                                             write_12st_fasta, write_3di_fasta,
                                             write_all_probs, write_fail_ids,
                                             write_mean_probs,
                                             write_profiles_text)
    from pholdlib.prostt5.device import parse_gpus

    from phold.databases.db import (check_modernprost_download,
                                    modernprost_zenodo_downloader)
    from phold.features.predict_3Di import write_embeddings

    spec = resolve_modernprost_model(model_name)
    task = str(task).lower()

    # ── resolve devices ─────────────────────────────────────────────────────
    devices = parse_gpus(cpu, gpus)
    logger.info(
        f"Beginning ModernProst ({spec.short_name}, task={task}) predictions "
        f"on device(s): {devices}"
    )
    if half_precision and devices == ["cpu"]:
        logger.info("CPU device — forcing full-precision (half-precision disabled).")
        half_precision = False
    logger.info(
        f"Using models in {'half' if half_precision else 'full'}-precision"
    )

    # ── flatten contigs into one globally sorted list ───────────────────────
    # Sorting globally means every multi-GPU shard starts length-descending,
    # so padding waste is balanced across devices. The nesting is restored
    # afterwards from cds_dict's own key order.
    flat_seq_dict: List[Tuple[str, str, int]] = []
    all_fail_ids: List[str] = []

    for record_id, seq_record_dict in cds_dict.items():
        for k, feat in seq_record_dict.items():
            v = feat.qualifiers.get("translation")
            if v and isinstance(v, str):
                flat_seq_dict.append((k, v, len(v)))
            else:
                logger.warning(
                    f"Protein header {k} is corrupt. It will be saved in fails.tsv"
                )
                all_fail_ids.append(k)

    flat_seq_dict.sort(key=lambda x: x[2], reverse=True)

    # ── inference ───────────────────────────────────────────────────────────
    result: ModernProstResult = run_modernprost_inference_multi_gpu(
        flat_seq_dict,
        devices=devices,
        model_dir=model_dir,
        model_name=spec.hf_name,
        half_precision=half_precision,
        threads=threads,
        task=task,
        check_fn=check_modernprost_download,
        zenodo_fn=modernprost_zenodo_downloader(spec.hf_name),
        max_residues=max_residues,
        max_seq_len=max_seq_len,
        max_batch=max_batch,
        chunk_len=DEFAULT_CHUNK_LEN if chunk_len is None else chunk_len,
        output_probs=True,  # phold always needs per-residue probs to mask AAs
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        desc="Predicting 3Di + 12st",
    )
    all_fail_ids.extend(result.fail_ids)

    # In the pssm task the engine returns distributions only. Collapse them to
    # the same (pred, mean_prob, all_prob) shape so the FASTA writers, the
    # confidence CSV and the amino-acid masking behave identically either way.
    if task == PSSM_TASK:
        flat_preds_3di = predictions_from_profiles(result.profiles_3di)
        flat_preds_12st = predictions_from_profiles(result.profiles_12st)
    else:
        flat_preds_3di = result.predictions_3di
        flat_preds_12st = result.predictions_12st

    # ── re-nest by contig, restoring the original per-contig key order ──────
    predictions_3di: Dict = {}
    predictions_12st: Dict = {}
    profiles: Dict = {}
    embeddings_per_residue: Dict = {}
    embeddings_per_protein: Dict = {}
    ordered_keys: List[str] = []
    header_map: Dict[str, str] = {}

    for record_id, seq_record_dict in cds_dict.items():
        keys = list(seq_record_dict.keys())
        predictions_3di[record_id] = {k: flat_preds_3di[k] for k in keys if k in flat_preds_3di}
        predictions_12st[record_id] = {
            k: flat_preds_12st[k] for k in keys if k in flat_preds_12st
        }
        if task == PSSM_TASK:
            profiles[record_id] = {
                k: {"3di": result.profiles_3di[k], "12st": result.profiles_12st[k]}
                for k in keys
                if k in result.profiles_3di and k in result.profiles_12st
            }
        if save_per_residue_embeddings:
            embeddings_per_residue[record_id] = {
                k: result.embeddings_per_residue[k]
                for k in keys
                if k in result.embeddings_per_residue
            }
        if save_per_protein_embeddings:
            embeddings_per_protein[record_id] = {
                k: result.embeddings_per_protein[k]
                for k in keys
                if k in result.embeddings_per_protein
            }

        for k in keys:
            ordered_keys.append(k)
            header_map[k] = k if proteins_flag else f"{record_id}:{k}"

    # ── write outputs ───────────────────────────────────────────────────────
    if all_fail_ids:
        write_fail_ids(all_fail_ids, Path(out_path) / "fails.tsv")

    # mask_threshold is deliberately not passed: an 'X' cannot be packed into
    # the combined 3Di+12st byte (see the module docstring).
    write_3di_fasta(flat_preds_3di, output_3di, ordered_keys, header_map=header_map)
    write_12st_fasta(flat_preds_12st, output_12st, ordered_keys, header_map=header_map)

    if save_per_residue_embeddings:
        write_embeddings(embeddings_per_residue, output_h5_per_residue)
    if save_per_protein_embeddings:
        write_embeddings(embeddings_per_protein, output_h5_per_protein)

    mean_probs_path = Path(out_path) / mean_probs_filename(prefix)
    with atomic_write_path(mean_probs_path) as tmp:
        write_mean_probs(flat_preds_3di, tmp, ordered_keys)

    mean_probs_12st_path = (
        Path(out_path) / f"{prefix}_{MODEL_PREFIX}_12st_mean_probabilities.csv"
    )
    with atomic_write_path(mean_probs_12st_path) as tmp:
        write_mean_probs(flat_preds_12st, tmp, ordered_keys)

    if output_probs:
        all_probs_path = (
            Path(out_path) / f"{prefix}_{MODEL_PREFIX}_3di_all_probabilities.json"
        )
        with atomic_write_path(all_probs_path) as tmp:
            write_all_probs(flat_preds_3di, tmp, ordered_keys)

    if task == PSSM_TASK:
        with atomic_write_path(Path(out_path) / f"{prefix}_profile_3di.txt") as tmp:
            write_profiles_text(result.profiles_3di, tmp, THREEDI_ALPHABET, ordered_keys)
        with atomic_write_path(Path(out_path) / f"{prefix}_profile_12st.txt") as tmp:
            write_profiles_text(result.profiles_12st, tmp, SS12_ALPHABET, ordered_keys)

    return predictions_3di, predictions_12st, profiles
