#!/usr/bin/env python3

import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import polars as pl
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from loguru import logger

from phold.features.create_foldseek_db import (
    generate_foldseek_db_from_aa_3di, generate_foldseek_db_from_structures)
from phold.features.run_foldseek import create_result_tsv, run_foldseek_search
from phold.io.handle_genbank import write_genbank
from phold.io.sub_db_outputs import create_sub_db_outputs
from phold.results.topfunction import (calculate_topfunctions_results,
                                       get_topcustom_hits,
                                       calculate_qcov_tcov,
                                       get_topfunctions)
from phold.utils.util import atomic_write_path, replace_pipe_in_fastq


def assign_annotation_confidence(row, structures: bool) -> str:
    """Classify a per-CDS row into a confidence tier.

    Pulled out of subcommand_compare so it's importable for unit tests and
    so the polars migration can rewrite this as a vectorised pl.when()
    chain without touching the surrounding orchestration.

    Tiers:
      - 'none'     : annotation_method == 'none'
      - 'pharokka' : annotation_method == 'pharokka'
      - 'high' / 'medium' / 'low' : per the coverage + fident + evalue
        (+ prostt5_confidence, when structures=False) heuristics described
        in the inline comment in subcommand_compare.
    """
    if row["annotation_method"] == "none":
        return "none"
    if row["annotation_method"] == "pharokka":
        return "pharokka"

    if structures:
        if (
            row["qCov"] > 0.8
            and row["tCov"] > 0.8
            and (row["fident"] > 0.3 or float(row["evalue"]) < 1e-10)
        ):
            return "high"
        if (row["qCov"] > 0.8 or row["tCov"] > 0.8) and (
            row["fident"] > 0.3 or float(row["evalue"]) < 1e-5
        ):
            return "medium"
        return "low"

    # not structures: prostt5_confidence factors in
    if (
        row["qCov"] > 0.8
        and row["tCov"] > 0.8
        and (
            row["fident"] > 0.3
            or row["prostt5_confidence"] > 60
            or float(row["evalue"]) < 1e-10
        )
    ):
        return "high"
    if (
        (row["qCov"] > 0.8 or row["tCov"] > 0.8)
        and (row["fident"] > 0.3 or 45 <= row["prostt5_confidence"] <= 60)
        and float(row["evalue"]) < 1e-5
    ):
        return "medium"
    return "low"


def assign_annotation_confidence_expr(structures: bool) -> pl.Expr:
    """Vectorised polars equivalent of ``assign_annotation_confidence``.

    Returns a polars expression that classifies every row of a DataFrame
    into a confidence tier. Equivalent to the scalar Python function on
    every row, but evaluated as a single ``pl.when().then()`` chain
    (no row-wise Python callback) — typically ~100× faster on large
    inputs and the byte-identical migration target.

    Args:
        structures: True for structure-based input (no prostt5_confidence
            available); False for the prostt5 path.

    Returns:
        A polars expression suitable for ``.with_columns(... .alias("annotation_confidence"))``.
    """
    evalue_f = pl.col("evalue").cast(pl.Float64)
    qCov = pl.col("qCov")
    tCov = pl.col("tCov")
    fident = pl.col("fident")
    method = pl.col("annotation_method")

    if structures:
        high = (
            (qCov > 0.8) & (tCov > 0.8)
            & ((fident > 0.3) | (evalue_f < 1e-10))
        )
        medium = (
            ((qCov > 0.8) | (tCov > 0.8))
            & ((fident > 0.3) | (evalue_f < 1e-5))
        )
    else:
        prostt5 = pl.col("prostt5_confidence")
        high = (
            (qCov > 0.8) & (tCov > 0.8)
            & ((fident > 0.3) | (prostt5 > 60) | (evalue_f < 1e-10))
        )
        medium = (
            ((qCov > 0.8) | (tCov > 0.8))
            & ((fident > 0.3) | ((prostt5 >= 45) & (prostt5 <= 60)))
            & (evalue_f < 1e-5)
        )

    return (
        pl.when(method == "none").then(pl.lit("none"))
        .when(method == "pharokka").then(pl.lit("pharokka"))
        .when(high).then(pl.lit("high"))
        .when(medium).then(pl.lit("medium"))
        .otherwise(pl.lit("low"))
    )


# ── per-contig CDS / function / sub-DB counts (writes _all_cds_functions.tsv) ─
# Output schema matches the original block in subcommand_compare: per
# contig, one "CDS" row, then 10 function-category rows, then 5 sub-DB
# rows — in that exact order, with 0-count rows emitted for categories
# that don't appear.

_FUNCTION_CATEGORIES: List[str] = [
    "connector",
    "DNA, RNA and nucleotide metabolism",
    "head and packaging",
    "integration and excision",
    "lysis",
    "moron, auxiliary metabolic gene and host takeover",
    "other",
    "tail",
    "transcription regulation",
    "unknown function",
]
_SUBDB_CATEGORIES: List[Tuple[str, str]] = [
    # (description_in_output, phrog_value_to_count)
    ("VFDB_Virulence_Factors", "vfdb"),
    ("CARD_AMR",                "card"),
    ("ACR_anti_crispr",         "acr"),
    ("Defensefinder",           "defensefinder"),
    ("Netflax",                 "netflax"),
]


def _write_function_counts_table(merged_df: pl.DataFrame, out_path: Path) -> None:
    """Build the per-contig function-count summary table and write as TSV.

    Output rows per contig (in this exact order):
      1. ``CDS`` row with total CDS count for the contig.
      2. 10 function rows (one per ``_FUNCTION_CATEGORIES`` value, in order;
         0 if absent).
      3. 5 sub-DB rows (one per ``_SUBDB_CATEGORIES`` entry, in order;
         counted on the ``phrog`` column).

    Implementation: single ``group_by("contig_id").agg(...)`` pass over
    ``merged_df`` computes CDS total + 15 boolean-sum aggregations at once.
    The previous per-contig ``filter`` then 16 ``.sum()`` calls scanned the
    frame ``len(contigs)`` times — replaced by one scan.
    """
    wide = merged_df.group_by("contig_id", maintain_order=True).agg([
        pl.len().alias("CDS"),
        *[(pl.col("function") == cat).sum().alias(cat) for cat in _FUNCTION_CATEGORIES],
        *[(pl.col("phrog") == phrog_value).sum().alias(description)
          for description, phrog_value in _SUBDB_CATEGORIES],
    ])

    rows: List[dict] = []
    for r in wide.iter_rows(named=True):
        contig = r["contig_id"]
        rows.append({"Description": "CDS", "Count": int(r["CDS"]), "Contig": contig})
        for category in _FUNCTION_CATEGORIES:
            rows.append({"Description": category, "Count": int(r[category]), "Contig": contig})
        for description, _ in _SUBDB_CATEGORIES:
            rows.append({"Description": description, "Count": int(r[description]), "Contig": contig})

    pl.DataFrame(rows).write_csv(out_path, separator="\t")


def subcommand_compare(
    gb_dict: Dict[str, Dict[str, Union[SeqRecord, SeqFeature]]],
    output: Path,
    threads: int,
    evalue: float,
    card_vfdb_evalue: float,
    sensitivity: float,
    database: Path,
    prefix: str,
    predictions_dir: Optional[Path],
    structures: bool,
    structure_dir: Optional[Path],
    logdir: Path,
    filter_structures: bool,
    remote_flag: bool,
    proteins_flag: bool,
    fasta_flag: bool,
    separate: bool,
    max_seqs: int,
    ultra_sensitive: bool,
    extra_foldseek_params: str,
    custom_db: str,
    foldseek_gpu: bool,
    restart: bool = False,
    clustered_db=False, # always False - keep the code for compatibility if I ever revert later, but clustered DBs were not better
    gpus: Optional[str] = None,
) -> bool:
    """
    Compare 3Di or PDB structures to the Phold DB

    Parameters:
        gb_dict (Dict[str, Dict[str, Union[SeqRecord, SeqFeature]]]): Nested dictionary containing genomic data.
        output (Path): Path to the output directory.
        threads (int): Number of threads to use.
        evalue (float): E-value threshold.
        card_vfdb_evalue (float): E-value threshold for CARD and VFDB databases.
        sensitivity (float): Sensitivity threshold.
        database (Path): Path to the reference database.
        prefix (str): Prefix for output files.
        predictions_dir (Optional[Path]): Path to the directory containing predictions.
        structures (bool): Flag indicating whether structures files are used.
        structure_dir (Optional[Path]): Path to the directory containing structures (.pdb or .cif) files.
        logdir (Path): Path to the directory for log files.
        filter_structuress (bool): Flag indicating whether to filter structure files.
        remote_flag (bool): Flag indicating whether the analysis is remote.
        proteins_flag (bool): Flag indicating whether proteins are used.
        fasta_flag (bool): Flag indicating whether FASTA files are used.
        separate (bool): Flag indicating whether to separate the analysis.
        max_seqs (int): Maximum results per query sequence allowed to pass the prefilter for foldseek.
        ultra_sensitive (bool): Whether to skip foldseek prefilter for maximum sensitivity
        extra_foldseek_params (str): Extra foldseek search parameters
        custom_db (str): Custom foldseek database
        foldseek_gpu (bool): Use Foldseek-GPU acceleration and ungappedprefilter
        restart (bool): Restart from foldseek_results.tsv

    Returns:
        bool: True if sub-databases are created successfully, False otherwise.
    """

    if predictions_dir is None and structures is False:
        logger.error(
            f"You did not specify --predictions_dir or --structures. Please check "
        )

    if structures and structure_dir is None:
        logger.error(
            f"You specified --structures but you did not specify --structure_dir. Please check "
        )

    if structure_dir and structures is False:
        logger.error(
            f"You specified --structure_dir but you did not specify --structures. Please check "
        )

    if proteins_flag is True:
        cds_dict = gb_dict
        non_cds_dict = {}
    else:
        # Create a nested dictionary to store CDS features by contig ID
        cds_dict = {}

        # makes the nested dictionary {contig_id:{cds_id: cds_feature}}
        for record_id, record in gb_dict.items():
            cds_dict[record_id] = {}

            for cds_feature in record.features:
                # for all pharokka genbanks, correct the bad phrog cats
                if fasta_flag is False:
                    if cds_feature.type == "CDS":
                        # first try Pharokka (ID) - do the updating
                        try:
                            # cds_id = cds_feature.qualifiers["ID"][0]
                            # update DNA, RNA and nucleotide metabolism from pharokka as it is broken as of 1.6.1
                            if "DNA" in cds_feature.qualifiers["function"][0]:
                                cds_feature.qualifiers["function"][
                                    0
                                ] = "DNA, RNA and nucleotide metabolism"
                                cds_feature.qualifiers["function"] = [
                                    cds_feature.qualifiers["function"][0]
                                ]  # Keep only the first element
                            # moron, auxiliary metabolic gene and host takeover as it is broken as of 1.6.1
                            if "moron" in cds_feature.qualifiers["function"][0]:
                                cds_feature.qualifiers["function"][
                                    0
                                ] = "moron, auxiliary metabolic gene and host takeover"
                                cds_feature.qualifiers["function"] = [
                                    cds_feature.qualifiers["function"][0]
                                ]  # Keep only the first element

                            # update No_PHROGs_HMM to No_PHROGs - if input is from pharokka --fast
                            if cds_feature.qualifiers["phrog"][0] == "No_PHROGs_HMM":
                                cds_feature.qualifiers["phrog"][0] = "No_PHROG"

                            cds_dict[record_id][
                                cds_feature.qualifiers["ID"][0]
                            ] = cds_feature

                        # not pharokka - must be from ncbi genbank (supported only)
                        except:
                            try:

                                # some NCBI Genbank CDS are actually pseudos
                                # e.g. OM418625

                                #  CDS             19638..19895
                                #                  /locus_tag="CPT_lambdaimm21_023"
                                #                  /pseudogene="unknown"
                                #                  /codon_start=1
                                #                  /transl_table=11
                                #                  /product="tail fiber protein stf"
                                if "pseudogene" in cds_feature.qualifiers:

                                    #logger.warning(f"Skipping pseudogene: {cds_feature}")
                                    continue
                                else:

                                    # add these extra fields to make it all play nice
                                    cds_feature.qualifiers["ID"] = cds_feature.qualifiers[
                                        "protein_id"
                                    ]
                                    cds_feature.qualifiers["function"] = []
                                    cds_feature.qualifiers["function"].append(
                                        "unknown function"
                                    )
                                    cds_feature.qualifiers["phrog"] = []
                                    cds_feature.qualifiers["phrog"].append("No_PHROG")

                                    cds_dict[record_id][
                                        cds_feature.qualifiers["ID"][0]
                                    ] = cds_feature

                            except:
                                logger.error(
                                    f"Feature {cds_feature} has no 'ID' or 'protein_id' qualifier in the Genbank file. Please add one in."
                                )

                # for fasta, will be fine to just add
                else:
                    cds_dict[record_id][cds_feature.qualifiers["ID"]] = cds_feature

        # non cds dict for later in genbank
        non_cds_dict = {}

        # makes the nested dictionary {contig_id:{non_cds_id: non_cds_feature}}
        for record_id, record in gb_dict.items():
            non_cds_dict[record_id] = {}

            i = 1
            for non_cds_feature in record.features:
                # Capture: anything not a CDS, OR a CDS marked as a pseudogene
                # (those don't have translations and must be handled as non-CDS).
                # NB: was previously `cds_feature.qualifiers` here — leftover
                # closure variable from the outer CDS loop. That checked the
                # wrong feature (the last CDS of the last record), causing
                # silent misclassification of pseudo-CDS, and a NameError on
                # records with zero CDS features.
                if non_cds_feature.type != "CDS" or (
                    non_cds_feature.type == "CDS"
                    and "pseudogene" in non_cds_feature.qualifiers
                ):
                    try:
                        non_cds_dict[record_id][
                            non_cds_feature.qualifiers["ID"][0]
                        ] = non_cds_feature
                    except:
                        non_cds_dict[record_id][
                            f"non_cds_feature_{i}"
                        ] = non_cds_feature
                        i += 1

    if not restart:

        # input predictions or structures
        if structures is False:
            # prostT5
            fasta_aa_input: Path = Path(predictions_dir) / f"{prefix}_aa.fasta"
            fasta_3di_input: Path = Path(predictions_dir) / f"{prefix}_3di.fasta"

        fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
        fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"

        ## copy the AA and 3Di from predictions directory if structures is false and phold compare is the command
        #
        # All three write paths below go via atomic_write_path so that an
        # interrupted predict-then-compare leaves the final files exactly
        # as they were (or absent) — never half-written. Previously the
        # raw ``shutil.copyfile`` + ``open(..., "w+")`` calls wrote
        # straight to the final paths; if the process died mid-copy or
        # mid-loop (OOM, Ctrl-C, disk full), the truncated AA / 3Di
        # FASTAs sat on disk and ``--restart`` happily picked them up as
        # if they were complete. atomic_write_path catches BaseException,
        # cleans up the sibling temp on failure, and only does the
        # ``os.replace`` swap onto the target on success.
        if structures is False:
            # if remote, these will not exist
            if remote_flag is False:
                if fasta_3di_input.exists():
                    logger.info(
                        f"Checked that the 3Di CDS file {fasta_3di_input} exists from phold predict"
                    )
                    with atomic_write_path(fasta_3di) as tmp:
                        shutil.copyfile(fasta_3di_input, tmp)
                else:
                    logger.error(
                        f"The 3Di CDS file {fasta_3di_input} does not exist. Please run phold predict and/or check the prediction directory {predictions_dir}"
                    )
                # copy the aa to file
                if fasta_aa_input.exists():
                    logger.info(
                        f"Checked that the AA CDS file {fasta_aa_input} exists from phold predict."
                    )
                    with atomic_write_path(fasta_aa) as tmp:
                        shutil.copyfile(fasta_aa_input, tmp)
                else:
                    logger.error(
                        f"The AA CDS file {fasta_aa_input} does not exist. Please run phold predict and/or check the prediction directory {predictions_dir}"
                    )
        ## write the AAs to file if structures is true because can't just copy from prediction_dir
        else:
            ## write the CDS to file
            logger.info(f"Writing the AAs to file {fasta_aa}.")
            # ``"w+"`` was unnecessary — the loop only writes — so this is
            # plain ``"w"`` on the sibling temp file.
            with atomic_write_path(fasta_aa) as tmp_fasta, open(tmp_fasta, "w") as out_f:
                for record_id, rest in cds_dict.items():
                    aa_contig_dict = cds_dict[record_id]

                    # writes the CDS to file
                    for seq_id, cds_feature in aa_contig_dict.items():
                        # if proteins, don't want the 'proteins:' as CDS id
                        if proteins_flag:
                            header = f">{seq_id}\n"
                            seq = f"{cds_feature.qualifiers['translation']}\n"
                        else:  # if genbank entry need to take the first seq as it is parsed as a list
                            header = f">{record_id}:{seq_id}\n"
                            seq = f"{cds_feature.qualifiers['translation'][0]}\n"
                        out_f.write(header)
                        out_f.write(seq)

        ############
        # create foldseek db
        ############

        foldseek_query_db_path: Path = Path(output) / "foldseek_db"
        foldseek_query_db_path.mkdir(parents=True, exist_ok=True)

        if structures is True:
            logger.info("Creating a foldseek query database from structures.")

            filtered_structures_path: Path = Path(output) / "filtered_structures"
            if filter_structures is True:
                logger.info(
                    f"--filter_structures is {filter_structures}. Therefore only the .pdb or .cif structure files with matching CDS ids will be copied and compared."
                )
                filtered_structures_path.mkdir(parents=True, exist_ok=True)

            generate_foldseek_db_from_structures(
                fasta_aa,
                foldseek_query_db_path,
                structure_dir,
                filtered_structures_path,
                logdir,
                prefix,
                filter_structures,
                proteins_flag,
            )
        else:
            generate_foldseek_db_from_aa_3di(
                fasta_aa, fasta_3di, foldseek_query_db_path, logdir, prefix
            )

        short_db_name = prefix

        # #  db search - not clustered

        database_name = "all_phold_structures"

        if clustered_db:
            database_name = "all_phold_structures_clustered_searchDB"



        if short_db_name == database_name:
            logger.error(
                f"Please choose a different {prefix} as this conflicts with the {database_name}"
            )

        #####
        # foldseek search
        #####


        query_db: Path = Path(foldseek_query_db_path) / short_db_name
        target_db: Path = Path(database) / database_name

        # make result and temp dirs
        result_db_base: Path = Path(output) / "result_db"
        result_db_base.mkdir(parents=True, exist_ok=True)
        result_db: Path = Path(result_db_base) / "result_db"

        temp_db: Path = Path(output) / "temp_db"
        temp_db.mkdir(parents=True, exist_ok=True)

        # make result tsv
        result_tsv: Path = Path(output) / "foldseek_results.tsv"

        # run foldseek search
        run_foldseek_search(
            query_db,
            target_db,
            result_db,
            temp_db,
            threads,
            logdir,
            evalue,
            sensitivity,
            max_seqs,
            ultra_sensitive,
            extra_foldseek_params,
            foldseek_gpu,
            structures,
            clustered_db,
            gpus=gpus,
        )

        
        create_result_tsv(query_db, target_db, result_db, result_tsv, logdir, foldseek_gpu, structures, threads)


   # restart

    ######
    # remove pipe in AA and 3Di FASTA (issue #86)
    # Won't exist if structures are used
    ######
    
    if not structures:
        fasta_aa: Path = Path(output) / f"{prefix}_aa.fasta"
        fasta_3di: Path = Path(output) / f"{prefix}_3di.fasta"

        replace_pipe_in_fastq(fasta_aa)
        replace_pipe_in_fastq(fasta_3di)


    ########
    # get topfunction
    ########

    # get top hit with non-unknown function for each CDS
    # also calculate the weighted bitscore df

    result_tsv: Path = Path(output) / "foldseek_results.tsv"
    database_name = "all_phold_structures"
    if clustered_db:
        database_name = "all_phold_structures_clustered_searchDB"

    filtered_topfunctions_df, weighted_bitscore_df = get_topfunctions(
        result_tsv, database, database_name, structures, card_vfdb_evalue, proteins_flag
    )

    # update the CDS dictionary with the tophits
    updated_cds_dict, filtered_tophits_df, source_dict = calculate_topfunctions_results(
        filtered_topfunctions_df,
        cds_dict,
        output,
        structures,
        proteins_flag,
        fasta_flag,
    )

    # generate per CDS foldseek information df and write to genbank
    per_cds_df = write_genbank(
        updated_cds_dict,
        non_cds_dict,
        source_dict,
        prefix,
        gb_dict,
        output,
        proteins_flag,
        separate,
        fasta_flag,
    )

    # ── normalise weighted_bitscore_df's query column → cds_id ────────────
    # For prostt5 input, ``query`` is "<contig_id>:<cds_id>" — split out
    # and keep cds_id. For structures / proteins input, query *is* the
    # cds_id, just rename.
    if not structures and not proteins_flag:
        weighted_bitscore_df = (
            weighted_bitscore_df
            .with_columns(
                pl.col("query")
                .str.splitn(":", 2)
                .struct.rename_fields(["_contig_id", "cds_id"])
                .alias("_q")
            )
            .unnest("_q")
            .drop(["query", "_contig_id"])
        )
    else:
        weighted_bitscore_df = weighted_bitscore_df.rename({"query": "cds_id"})

    # Drop columns that would otherwise duplicate with per_cds_df during
    # the subsequent left-join. ``query`` is dropped because cds_id is the
    # downstream identifier — keeping ``query`` would add a redundant col
    # (e.g. "<contig>:<cds_id>") to per_cds_predictions.tsv.
    drop_dupes = ["query", "phrog", "product", "function"]
    if not proteins_flag:
        drop_dupes = ["contig_id"] + drop_dupes
    drop_dupes = [c for c in drop_dupes if c in filtered_tophits_df.columns]
    filtered_tophits_df = filtered_tophits_df.drop(drop_dupes)

    # Two left-joins to merge tophits + weighted bitscores onto per_cds_df.
    merged_df = (
        per_cds_df
        .join(filtered_tophits_df,    on="cds_id", how="left")
        .join(weighted_bitscore_df,   on="cds_id", how="left")
    )

    # Pandas-parity float casts. The original pandas implementation
    # left-joined int-typed foldseek/weighted columns onto rows that may
    # not have a match; pandas upcasts the column to float64 to accommodate
    # the NaN padding, so the TSV ends up with "2019.0" / "1.0" / "0.0".
    # Polars supports nullable ints natively and keeps the column as Int64
    # (writes "2019" / "1" / "0"). Cast explicitly so per_cds_predictions.tsv
    # stays byte-identical to the pre-migration output (and the CLI
    # comparison suite keeps passing on exact diff).
    _pandas_parity_float_cols = [
        # from filtered_tophits_df
        "bitscore", "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        # from weighted_bitscore_df (the bitscore_proportion columns that
        # stayed Int64 in topfunction.py's pandas-quirk preservation;
        # they'd be Float64 after the pandas left-join upcast).
        "integration_and_excision_bitscore_proportion",
        "connector_bitscore_proportion",
        "lysis_bitscore_proportion",
        "other_bitscore_proportion",
        "unknown_function_bitscore_proportion",
        "moron_auxiliary_metabolic_gene_and_host_takeover_bitscore_proportion",
        "transcription_regulation_bitscore_proportion",
    ]
    merged_df = merged_df.with_columns(
        [pl.col(c).cast(pl.Float64) for c in _pandas_parity_float_cols if c in merged_df.columns]
    )

    # Move ``annotation_method`` to immediately after ``product``. The
    # proteins-compare path produces a merged_df without that column —
    # pandas' .reindex(columns=) used to silently fill it with NaN, so
    # we replicate that by adding a null column when missing (polars'
    # .select() would otherwise raise ColumnNotFoundError).
    if "annotation_method" not in merged_df.columns:
        merged_df = merged_df.with_columns(
            pl.lit(None, dtype=pl.Utf8).alias("annotation_method")
        )
    cols = merged_df.columns
    product_idx = cols.index("product")
    new_order = (
        [c for c in cols[: product_idx + 1] if c != "annotation_method"]
        + ["annotation_method"]
        + [c for c in cols[product_idx + 1 :] if c != "annotation_method"]
    )
    merged_df = merged_df.select(new_order)

    # add qcov and tcov
    merged_df = calculate_qcov_tcov(merged_df)

    # NEEDS TO READ IN THE f"{prefix}_prostT5_3di_mean_probabilities.csv" - can't pass from the predict function in case using phold compare
    if predictions_dir is None:  # if running phold run
        mean_probs_out_path: Path = (
            Path(output) / f"{prefix}_prostT5_3di_mean_probabilities.csv"
        )
    else:  # if running phold compare or phold proteins-compare
        mean_probs_out_path: Path = (
            Path(predictions_dir) / f"{prefix}_prostT5_3di_mean_probabilities.csv"
        )

    # Merge in ProstT5 confidence scores — only for the non-structures path.
    if not structures:
        prostT5_conf_df = pl.read_csv(
            mean_probs_out_path,
            separator=",",
            has_header=False,
            new_columns=["cds_id", "prostt5_confidence"],
        )
        merged_df = merged_df.join(prostT5_conf_df, on="cds_id", how="left")

    # ── confidence classification ──────────────────────────────────────────
    # High - 80%+ reciprocal coverage + one of i) >30% seqid cutoff for the
    #        light zone (https://doi.org/10.1093/protein/12.2.85) OR
    #        ii) ProstT5 confidence > 60% OR evalue < 1e-10
    # Medium - either qCov or tCov 80%+ + (30%+ seqid OR ProstT5 confidence
    #          in [45,60]) AND evalue < 1e-5
    # Low - everything else
    # Vectorised pl.when() chain replaces the row-wise apply.
    merged_df = merged_df.with_columns(
        assign_annotation_confidence_expr(structures=structures)
        .alias("annotation_confidence")
    )

    # Move ``annotation_confidence`` to immediately after ``annotation_method``.
    cols = merged_df.columns
    method_idx = cols.index("annotation_method")
    new_order = (
        [c for c in cols[: method_idx + 1] if c != "annotation_confidence"]
        + ["annotation_confidence"]
        + [c for c in cols[method_idx + 1 :] if c != "annotation_confidence"]
    )
    merged_df = merged_df.select(new_order)

    # save
    merged_df_path: Path = Path(output) / f"{prefix}_per_cds_predictions.tsv"
    merged_df.write_csv(merged_df_path, separator="\t")

    # custom db output

    #####
    # custom db
    #####

    if custom_db:

        logger.info(
            f"Foldseek will also be run against your custom database {custom_db}"
        )
        # make result and temp dirs
        result_db_custom: Path = Path(result_db_base) / "result_db_custom"
        result_tsv_custom: Path = Path(output) / "foldseek_results_custom.tsv"

        # try:
        run_foldseek_search(
            query_db,
            Path(custom_db),
            result_db_custom,
            temp_db,
            threads,
            logdir,
            evalue,
            sensitivity,
            max_seqs,
            ultra_sensitive,
            extra_foldseek_params,
            foldseek_gpu,
            structures,
            clustered_db=False,  # no custom db cluster searching
            gpus=gpus,
        )

        # make result tsv
        create_result_tsv(
            query_db,
            Path(custom_db),
            result_db_custom,
            result_tsv_custom,
            logdir,
            foldseek_gpu,
            structures,
            threads,
        )

        tophit_custom_df = get_topcustom_hits(
            result_tsv_custom, structures, proteins_flag
        )

        #### merge
        # left merge on the cds_id to get every cds/contig id (make it easier for downstream processing)

        # Left-join custom hits onto every CDS so downstream code sees
        # every cds_id (possibly with NULL custom-hit fields). cds_id is
        # always unique.
        all_cds_cols = ["cds_id"] if proteins_flag else ["contig_id", "cds_id"]
        all_cds_df = merged_df.select(all_cds_cols)
        tophit_custom_df = all_cds_df.join(tophit_custom_df, how="left", on="cds_id")

        # Reorder columns: cds_id (and contig_id for non-proteins) first.
        front_cols = ["cds_id"] if proteins_flag else ["contig_id", "cds_id"]
        rest = [c for c in tophit_custom_df.columns if c not in ("contig_id", "cds_id")]
        tophit_custom_df = tophit_custom_df.select(front_cols + rest)

        # get coverages
        tophit_custom_df = calculate_qcov_tcov(tophit_custom_df)
        custom_hits_path: Path = Path(output) / f"{prefix}_custom_database_hits.tsv"
        tophit_custom_df.write_csv(custom_hits_path, separator="\t")

        # except:
        #     logger.error(f"Foldseek failed to run against your custom database {custom_db}. Please check that it is formatted correctly as a Foldseek database")

    # sub dbs output
    # save vfdb card acr defensefinder hits with more metadata
    sub_dbs_created = create_sub_db_outputs(merged_df, database, output)

    # Per-contig CDS / function / sub-DB count summary written to
    # <prefix>_all_cds_functions.tsv. Replaces the original ~100-line
    # block of repeated pd.DataFrame constructions + concat.
    if not proteins_flag:
        descriptions_total_path: Path = Path(output) / f"{prefix}_all_cds_functions.tsv"
        _write_function_counts_table(merged_df, descriptions_total_path)

    return sub_dbs_created
