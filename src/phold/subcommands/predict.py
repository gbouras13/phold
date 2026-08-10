#!/usr/bin/env python3

from pathlib import Path
from typing import Optional

import numpy as np
from loguru import logger

from pholdlib.databases.modernprost import (
    CLASSIFICATION_TASK,
    PSSM_TASK,
    default_task_for,
    is_modernprost_model,
)

from phold.features.predict_3Di import get_embeddings
from phold.features.predict_3di_12st import get_modernprost_predictions


def mask_low_confidence_aa(sequence: str, scores, threshold: float = 0.5) -> str:
    """
    Replace amino acids whose ProstT5 confidence score is below *threshold* with 'X'.

    *scores* is a numpy array of shape (1, L) or (L,), or an equivalent nested list.
    """
    score_arr = np.asarray(scores, dtype=np.float64).ravel()
    mask = score_arr < threshold
    if not mask.any():
        return sequence
    buf = np.frombuffer(bytearray(sequence, "ascii"), dtype=np.uint8).copy()
    buf[mask] = 88  # ord("X")
    return buf.tobytes().decode("ascii")


def subcommand_predict(
    gb_dict: dict,
    method: str,
    output: Path,
    prefix: str,
    cpu: bool,
    omit_probs: bool,
    model_dir: Path,
    model_name: str,
    checkpoint_path: Path,
    batch_size: int,
    proteins_flag: bool,
    fasta_flag: bool,
    save_per_residue_embeddings: bool,
    save_per_protein_embeddings: bool,
    threads: int,
    mask_threshold: float,
    hyps: bool,
    gpus: Optional[str] = None,
    task: str = "auto",
    max_batch_residues: Optional[int] = None,
    logdir: Optional[Path] = None,
) -> bool:
    """
    Wrapper command for phold predict.

    Runs either ProstT5 (encoder + CNN prediction head, 3Di only) or a
    ModernProst checkpoint (single encoder, 3Di + 12-state), selected by
    *model_name*.

    Args:
        gb_dict (Dict[str, any]): Dictionary containing GenBank records.
        method (str): "pharokka", "bakta" or "ncbi" - input format
        output (str): Output directory path.
        prefix (str): Prefix for output file names.
        cpu (bool): Flag indicating whether to use CPU for prediction.
        omit_probs (bool): Flag indicating whether to omit per-residue prediction probabilities.
        model_dir (str): Directory containing the model.
        model_name (str): HuggingFace id or ModernProst registry name.
        checkpoint_path (Path): Path to the ProstT5 CNN checkpoint. Unused for ModernProst.
        batch_size (int): Batch size for prediction.
        proteins_flag (bool): True if phold proteins-predict, false otherwise
        fasta_flag (bool): True if pyrodigal-gv was used to predict CDS from FASTA input. False otherwise
        save_per_residue_embeddings (bool, optional): Whether to save per residue embeddings to h5 file. Defaults to False.
        save_per_protein_embeddings (bool, optional): Whether to save mean per protein embeddings to h5 file. Defaults to False.
        threads (int): number of cpu threads
        task (str): "auto", "classification" or "pssm". ModernProst only; "auto"
            uses the task the checkpoint was trained for.
        max_batch_residues (Optional[int]): Max residues per inference batch.
            Defaults to the per-backend value when None.
        logdir (Optional[Path]): Log directory, needed to build the Foldseek
            profile database in the pssm task.

    Returns:
        bool: True if prediction succeeds, False otherwise.
    """

    #########
    # make nested dictionary
    #########

    if hyps:
        if method == "Pharokka":
            logger.info(
                f"You have used --hyps and a Pharokka style input Genbank was detected."
            )
            logger.info(
                "Only unknown function proteins from your Pharokka input Genbank will be extracted and annotated with Phold."
            )
        else:
            logger.warning(
                "You can specified --hyps but your input Genbank file is not a Pharokka style input Genbank file."
            )
            logger.warning(
                "Ignoring --hyps: all input CDS will be annotated with Phold."
            )

    fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"


    # if proteins, already done and passed as gb_dict
    if proteins_flag is True:
        cds_dict = gb_dict
    else:
        # Create a nested dictionary to store CDS features by contig ID
        cds_dict = {}
        # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
        for record_id, record in gb_dict.items():
            cds_dict[record_id] = {}
            for cds_feature in record.features:
                if cds_feature.type == "CDS":
                    # due to the weird list issue when parsing from genbank file
                    if fasta_flag is False:

                        # cds_feature.qualifiers["translation"] = cds_feature.qualifiers[
                        #     "translation"
                        # ][0]

                        # some NCBI Genbank CDS are actually pseudos
                        # e.g. OM418625

                        #  CDS             19638..19895
                        #                  /locus_tag="CPT_lambdaimm21_023"
                        #                  /pseudogene="unknown"
                        #                  /codon_start=1
                        #                  /transl_table=11
                        #                  /product="tail fiber protein stf"

                        if (
                            "translation" not in cds_feature.qualifiers
                            or len(cds_feature.qualifiers["translation"]) == 0
                        ):
                            logger.warning(
                                f"Skipping CDS without provided translation in input, likely a pseudogene"
                            )
                            logger.warning(f"CDS: {cds_feature}"
                            )
                            continue

                        cds_feature.qualifiers["translation"] = cds_feature.qualifiers[
                            "translation"
                        ][0]

                        if method == "Pharokka":
                            try:
                                cds_id = cds_feature.qualifiers["ID"][
                                    0
                                ]  # if this breaks, will mean not Pharokka input

                                if hyps:
                                    if (
                                        cds_feature.qualifiers["function"][0]
                                        != "unknown function"
                                    ):
                                        logger.info(
                                            f"Skipping {cds_id} as it has a known function from Pharokka"
                                        )
                                        continue
                            except:
                                logger.error(
                                    f"Feature {cds_feature} has no 'ID' qualifier in the Genbank file despite being likely Pharokka origin. Please check your input Genbank file."
                                )
                        else:
                            # next try Genbank/NCBI (uses protein_id)
                            if method == "NCBI":
                                try:
                                    # add these extra fields to make it all play nice
                                    cds_feature.qualifiers["ID"] = (
                                        cds_feature.qualifiers["protein_id"]
                                    )
                                    cds_feature.qualifiers["function"] = []
                                    cds_feature.qualifiers["function"].append(
                                        "unknown function"
                                    )
                                    cds_feature.qualifiers["phrog"] = []
                                    cds_feature.qualifiers["phrog"].append("No_PHROG")

                                    cds_id = cds_feature.qualifiers["ID"][0]

                                    cds_dict[record_id][
                                        cds_feature.qualifiers["ID"][0]
                                    ] = cds_feature
                                except:
                                    logger.error(
                                        f"Feature {cds_feature} has no 'protein_ID' qualifier in the Genbank file despite being detected as being likely NCBI Refseq style. Please add one in."
                                    )
                            # finally try bakta (use locus_tag)
                            if method == "Bakta":
                                try:
                                    # add these extra fields to make it all play nice
                                    cds_feature.qualifiers["ID"] = (
                                        cds_feature.qualifiers["locus_tag"]
                                    )
                                    cds_feature.qualifiers["function"] = []
                                    cds_feature.qualifiers["function"].append(
                                        "unknown function"
                                    )
                                    cds_feature.qualifiers["phrog"] = []
                                    cds_feature.qualifiers["phrog"].append("No_PHROG")

                                    cds_id = cds_feature.qualifiers["locus_tag"][0]

                                    cds_dict[record_id][
                                        cds_feature.qualifiers["locus_tag"][0]
                                    ] = cds_feature
                                except:
                                    logger.error(
                                        f"Feature {cds_feature} has no 'locus_tag' qualifier in the Genbank file despite being detected as being likely bakta origin. Please check your input Genbank file."
                                    )
                        # append CDS
                        cds_dict[record_id][cds_id] = cds_feature

                    else:
                        cds_dict[record_id][cds_feature.qualifiers["ID"]] = cds_feature

    # issue #86 with GenBank format

    new_cds_dict = {}
    for record_id, record in cds_dict.items():
        if "~PIPE~" in record_id:
            logger.error(
                f"Your FASTA header {record_id} has __PIPE__ in the header. "
                "Please remove all instances of __PIPE__ in the header before running Phold "
                "(or Pharokka before Phold)"
            )
        else:
            record_id = record_id.replace("|", "~PIPE~")
        new_cds_dict[record_id] = record

    cds_dict = new_cds_dict
    del new_cds_dict


    fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"
    fasta_12st: Path = Path(output) / f"{prefix}_12st.fasta"
    # embeddings h5 - will only be generated if flag is true
    output_h5_per_residue: Path = Path(output) / f"{prefix}_embeddings_per_residue.h5"
    output_h5_per_protein: Path = Path(output) / f"{prefix}_embeddings_per_protein.h5"

    if cpu is True:
        half_precision = False
    else:
        half_precision = True

    if omit_probs:
        output_probs = False
    else:
        output_probs = True

    modernprost = is_modernprost_model(model_name)
    profiles = {}

    if modernprost:
        resolved_task = (
            default_task_for(model_name) if task in (None, "auto") else str(task).lower()
        )
        if resolved_task not in (CLASSIFICATION_TASK, PSSM_TASK):
            logger.error(
                f"--task must be auto, {CLASSIFICATION_TASK} or {PSSM_TASK}; "
                f"got {task!r}"
            )
        if resolved_task != default_task_for(model_name):
            logger.warning(
                f"{model_name} was trained for the "
                f"'{default_task_for(model_name)}' task but --task "
                f"{resolved_task} was requested."
            )

        if mask_threshold and mask_threshold > 0:
            # The combined 3Di+12st Foldseek encoding has no masked state, so
            # the 3Di FASTA is written unmasked. Masking still applies to the
            # amino-acid FASTA below.
            logger.info(
                f"--mask_threshold {mask_threshold} applies to the amino-acid "
                "FASTA only when using a ModernProst model: the combined "
                "3Di+12st Foldseek alphabet cannot represent a masked residue."
            )

        predictions, predictions_12st, profiles = get_modernprost_predictions(
            cds_dict,
            output,
            prefix,
            model_dir,
            model_name,
            fasta_3di,
            fasta_12st,
            output_h5_per_residue,
            output_h5_per_protein,
            half_precision=half_precision,
            task=resolved_task,
            max_residues=max_batch_residues if max_batch_residues else 50000,
            max_seq_len=30000,
            max_batch=batch_size,
            cpu=cpu,
            output_probs=output_probs,
            proteins_flag=proteins_flag,
            save_per_residue_embeddings=save_per_residue_embeddings,
            save_per_protein_embeddings=save_per_protein_embeddings,
            threads=threads,
            gpus=gpus,
        )
    else:
        if task not in (None, "auto"):
            logger.warning(
                f"--task {task} only applies to the ModernProst models; "
                f"ignoring it for {model_name}."
            )
        predictions = get_embeddings(
            cds_dict,
            output,
            prefix,
            model_dir,
            model_name,
            checkpoint_path,
            fasta_3di,
            output_h5_per_residue,
            output_h5_per_protein,
            half_precision=half_precision,
            max_residues=max_batch_residues if max_batch_residues else 5000,
            max_seq_len=1000,
            max_batch=batch_size,
            cpu=cpu,
            output_probs=output_probs,
            proteins_flag=proteins_flag,
            save_per_residue_embeddings=save_per_residue_embeddings,
            save_per_protein_embeddings=save_per_protein_embeddings,
            threads=threads,
            mask_threshold=mask_threshold,
            gpus=gpus,
        )

    mask_prop_threshold = mask_threshold / 100

    ########
    ## write the AA CDS to file
    ######

    with open(fasta_aa, "w") as out_f:
        for contig_id, rest in cds_dict.items():
            aa_contig_dict = cds_dict[contig_id]
            prediction_contig_dict = predictions[contig_id]
            prediction_contig_dict = {
                k: v for k, v in prediction_contig_dict.items() if len(v[0]) > 0
            }
            parts = []
            for seq_id, cds_feature in aa_contig_dict.items():
                header = f">{seq_id}\n" if proteins_flag else f">{contig_id}:{seq_id}\n"

                prot_seq = cds_feature.qualifiers["translation"]

                try:
                    # this will fail if ProstT5 OOM fails (or fails for some other reason)
                    prot_seq = mask_low_confidence_aa(
                        prot_seq,
                        prediction_contig_dict[seq_id][2],
                        threshold=mask_prop_threshold,
                    )
                except (KeyError, IndexError):
                    # in that case, just return 'X' aka masked proteins
                    prot_seq = "X" * len(prot_seq)

                parts.append(f"{header}{prot_seq}\n")
            out_f.write("".join(parts))

    ########
    ## build the Foldseek query profile DB (pssm task only)
    ######
    # Deliberately after the amino-acid FASTA is written: the profile DB reads
    # its sequences from that file, so the amino-acid channel of the profile
    # matches the masking that was actually applied.
    if profiles:
        from phold.features.create_foldseek_db import generate_foldseek_profile_db

        profile_db_path: Path = Path(output) / "query_profiledb"
        logger.info(f"Building Foldseek query profile databases in {profile_db_path}")
        generate_foldseek_profile_db(
            profiles,
            fasta_aa,
            profile_db_path,
            logdir if logdir is not None else Path(output) / "logs",
            prefix,
        )

    return True
