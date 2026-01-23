#!/usr/bin/env python3

from pathlib import Path

from loguru import logger

from phold.features.predict_3Di import get_embeddings
from phold.features.profiles import generate_mmseqs_db_from_aa, build_lookup

def mask_low_confidence_aa(sequence, scores, threshold=0.25):
    """
    Masks all low confidence AA to X if their corresponding ProstT5 confidence score is below the given threshold.

    Parameters:
    sequence (str): The amino acid sequence.
    scores (List[float]): A list of confidence scores for each amino acid.
    threshold (float, optional): The confidence threshold below which amino acids are converted to lowercase. Default is 0.5.

    Returns:
    str: The modified amino acid sequence with low-confidence residues in lowercase.
    """
    return "".join(
        "X" if float(score) < threshold else aa # scores is a list so is fine no need to unpack with *scores
        for aa, score in zip(sequence, scores)
    )


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
    max_batch_residues: int,
    batch_size: int,
    proteins_flag: bool,
    fasta_flag: bool,
    save_per_residue_embeddings: bool,
    save_per_protein_embeddings: bool,
    threads: int,
    mask_threshold: float,
    hyps: bool
) -> bool:
    """
    Wrapper command for phold predict. Predicts embeddings using ProstT5 encoder + CNN prediction head.

    Args:
        gb_dict (Dict[str, any]): Dictionary containing GenBank records.
        method (str): "pharokka", "bakta" or "ncbi" - input format
        output (str): Output directory path.
        prefix (str): Prefix for output file names.
        cpu (bool): Flag indicating whether to use CPU for prediction.
        omit_probs (bool): Flag indicating whether to omit prediction probabilities from ProstT5.
        model_dir (str): Directory containing the ProstT5 model.
        model_name (str): Name of the ProstT5 model.
        checkpoint_path (Path): Path to ProstT5 CNN checkpoint.
        max_batch_residues (int): Maximum total residues per batch
        batch_size (int): Batch size for prediction.
        proteins_flag (bool): True if phold proteins-predict, false otherwise
        fasta_flag (bool): True if pyrodigal-gv was used to predict CDS from FASTA input. False otherwise
        save_per_residue_embeddings (bool, optional): Whether to save per residue embeddings to h5 file. Defaults to False.
        save_per_protein_embeddings (bool, optional): Whether to save mean per protein embeddings to h5 file. Defaults to False.
        threads (int): number of cpu threads

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

    # for profiles - build the MMseqs2 DB
    logdir: Path = Path(output) / "logs"

    if model_name == "gbouras13/modernprost-profiles":


        lookup = generate_mmseqs_db_from_aa(cds_dict, output, logdir, prefix, proteins_flag)

        mmseqs2_db_dir: Path = Path(output) / f"query_profiledb" # this will be subdir where the dummy mmseqs2 db is 
        mmseqs2_db_dir.mkdir(parents=True, exist_ok=True)

        # create MMSeqs2 db names
        short_db_name = f"{prefix}"

        mmseqs_db: Path = Path(mmseqs2_db_dir) / f"{short_db_name}"

        mmseqs_db_lookup: Path = Path(mmseqs2_db_dir) / f"{short_db_name}.lookup"

        seq_output_db: Path = Path(mmseqs2_db_dir) / f"{short_db_name}_profile"
        seq_index = f"{seq_output_db}.index"
        pssm_output_db: Path = Path(mmseqs2_db_dir) / f"{short_db_name}_profile_ss"
        pssm_index = f"{pssm_output_db}.index"
        
        # lookup = build_lookup(mmseqs_db_lookup)

        pssm_db_builder = {
            "dummy_seq_db": mmseqs_db,
            "lookup": lookup,
            "parts": [],
            "index_lines": [],
            "offset": 0,
            "output_db": pssm_output_db,
            "index_file": pssm_index,
            "seq_output_db": seq_output_db,
            "seq_index_file": seq_index
        }
    else:
        pssm_db_builder = None # for the function

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
        max_batch_residues=max_batch_residues,
        max_seq_len=10000,
        max_batch=batch_size,
        cpu=cpu,
        output_probs=output_probs,
        proteins_flag=proteins_flag,
        save_per_residue_embeddings=save_per_residue_embeddings,
        save_per_protein_embeddings=save_per_protein_embeddings,
        threads=threads,
        mask_threshold=mask_threshold,
        pssm_db_builder=pssm_db_builder
    )

    mask_prop_threshold = mask_threshold / 100

    if model_name == "gbouras13/modernprost-profiles":
        mask_prop_threshold = 0 # no masking for profiles

    ########
    ## write the AA CDS to file
    ######

    with open(fasta_aa, "w+") as out_f:
        for contig_id, rest in cds_dict.items():
            aa_contig_dict = cds_dict[contig_id]
            prediction_contig_dict = predictions[contig_id]
            prediction_contig_dict = {
                k: v for k, v in prediction_contig_dict.items() if len(v[0]) > 0
            }
            for seq_id, cds_feature in aa_contig_dict.items():
                if proteins_flag is True:
                    out_f.write(f">{seq_id}\n")
                else:
                    out_f.write(f">{contig_id}:{seq_id}\n")

                prot_seq = cds_feature.qualifiers["translation"]
                # prediction_contig_dict[seq_id][2] these are teh ProstT5 confidence scores from 0-1 - need to convert to list

                try:
                    # this will fail if ProstT5 OOM fails (or fails for some other reason)
                    prot_seq = mask_low_confidence_aa(
                        prot_seq,
                        prediction_contig_dict[seq_id][2].tolist(),
                        threshold=mask_prop_threshold,
                    )
                except (KeyError, IndexError):
                    # in that case, just return 'X' aka masked proteins
                    prot_seq = "X" * len(prot_seq)

                out_f.write(f"{prot_seq}\n")

    return True
