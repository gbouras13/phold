#!/usr/bin/env python3
from pathlib import Path
from typing import Dict, List, Tuple, Union

import polars as pl
from loguru import logger


# ── weighted-bitscore output schema ─────────────────────────────────────────
# (function-name-in-data, output-column-name) pairs. Order matters: it
# determines the column order in the weighted_bitscore output TSV.
_BITSCORE_CATEGORIES: List[Tuple[str, str]] = [
    ("head and packaging",
        "head_and_packaging_bitscore_proportion"),
    ("integration and excision",
        "integration_and_excision_bitscore_proportion"),
    ("tail",
        "tail_bitscore_proportion"),
    ("moron, auxiliary metabolic gene and host takeover",
        "moron_auxiliary_metabolic_gene_and_host_takeover_bitscore_proportion"),
    ("DNA, RNA and nucleotide metabolism",
        "DNA_RNA_and_nucleotide_metabolism_bitscore_proportion"),
    ("connector",
        "connector_bitscore_proportion"),
    ("transcription regulation",
        "transcription_regulation_bitscore_proportion"),
    ("lysis",
        "lysis_bitscore_proportion"),
    ("other",
        "other_bitscore_proportion"),
    # ``unknown function`` is always 0 in the output (the original loop skips
    # it when building weighted_counts_normalised, and .get default = 0).
    ("unknown function",
        "unknown_function_bitscore_proportion"),
]


def get_topfunctions(
    result_tsv: Path,
    database: Path,
    database_name: str,
    structures: bool,
    card_vfdb_evalue: float,
    proteins_flag: bool,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Process Foldseek output to extract top functions and weighted bitscores.

    Args:
        result_tsv (Path): Path to the Foldseek result TSV file.
        database (Path): Path to the database directory.
        database_name (str): Name of the database.
        structures (bool): Flag indicating whether structures have been added.
        card_vfdb_evalue (float): E-value threshold for card and vfdb hits.
        proteins_flag (bool): Flag indicating whether proteins are used.

    Returns:
        Tuple[pl.DataFrame, pl.DataFrame]: A tuple containing two DataFrames:
            1. ``topfunction_df`` — one best hit per query.
            2. ``weighted_bitscore_df`` — per-query proportion of bitscore
               accounted for by each functional category.
    """
    logger.info("Processing Foldseek output")

    base_cols = [
        "query", "target", "bitscore", "fident", "evalue",
        "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
    ]
    col_list = base_cols + (["alntmscore", "lddt"] if structures else [])

    # Force evalue to Utf8 — pandas had inferred it as object/string; the
    # ``"{:.3e}".format(float(x))`` formatting step downstream relies on
    # string input and the existing snapshots are pinned to string-style
    # scientific repr.
    foldseek_df = pl.read_csv(
        result_tsv,
        separator="\t",
        has_header=False,
        new_columns=col_list,
        schema_overrides={"evalue": pl.Utf8},
        infer_schema_length=10_000,
    )

    # in case the foldseek output is empty
    if foldseek_df.is_empty():
        logger.error(
            "Foldseek found no hits whatsoever - please check whether your "
            "input is really phage-like"
        )

    # issue #86 — restore the literal '|' in queries (mangled to '~PIPE~'
    # upstream to survive foldseek's tab-separated parsing).
    foldseek_df = foldseek_df.with_columns(
        pl.col("query").str.replace_all("~PIPE~", "|", literal=True),
        pl.lit("foldseek").alias("annotation_source"),
    )

    # Split out contig_id + cds_id (or just cds_id) from the query column.
    if not structures and not proteins_flag:
        # prostt5 path: query = "<contig_id>:<cds_id>" — split on the first ":"
        foldseek_df = foldseek_df.with_columns(
            pl.col("query")
            .str.splitn(":", 2)
            .struct.rename_fields(["contig_id", "cds_id"])
            .alias("_q")
        ).unnest("_q")
    else:
        # structures or proteins_flag path: query is the cds_id (possibly
        # with a ``.pdb`` / ``.cif`` suffix). NOTE — the original pandas
        # code re-assigns cds_id twice: the second assignment OVERWRITES
        # the first, replacing only ``.cif`` from the *original* query
        # (the ``.pdb`` replacement is silently discarded). Preserve that
        # behaviour for byte-identity with existing CLI outputs.
        foldseek_df = foldseek_df.with_columns(
            pl.col("query").str.replace_all(".cif", "", literal=True).alias("cds_id")
        )

    # Clean up ``.pdb`` suffix from target, then split target into phrog +
    # tophit_protein. Note: matches pandas which only strips ``.pdb`` (not
    # ``.cif``) from target. Preserved for byte-identity.
    foldseek_df = foldseek_df.with_columns(
        pl.col("target").str.replace_all(".pdb", "", literal=True)
    ).with_columns(
        pl.col("target")
        .str.splitn(":", 2)
        .struct.rename_fields(["phrog", "tophit_protein"])
        .alias("_t")
    ).unnest("_t").drop("target")

    # Strip the leading "phrog_" prefix from all phrog values.
    foldseek_df = foldseek_df.with_columns(
        pl.col("phrog").str.replace_all("phrog_", "", literal=True)
    )

    # envhog_ rows: strip the prefix from phrog and prepend it to tophit_protein.
    # Both expressions in this with_columns() see the SAME input — the
    # ``starts_with`` mask is evaluated against the unmodified phrog column
    # for both branches, so the rename and the prepend stay aligned.
    foldseek_df = foldseek_df.with_columns(
        pl.when(pl.col("phrog").str.starts_with("envhog_"))
          .then(pl.col("phrog").str.replace_all("envhog_", "", literal=True))
          .otherwise(pl.col("phrog"))
          .alias("phrog"),
        pl.when(pl.col("phrog").str.starts_with("envhog_"))
          .then(pl.lit("envhog_") + pl.col("tophit_protein"))
          .otherwise(pl.col("tophit_protein"))
          .alias("tophit_protein"),
    )

    # efam_ and dgr_ rows: just strip the prefix from phrog (no protein
    # rename). These run sequentially because each sees the phrog column
    # after the previous strip.
    foldseek_df = foldseek_df.with_columns(
        pl.col("phrog").str.replace_all("efam_", "", literal=True)
    ).with_columns(
        pl.col("phrog").str.replace_all("dgr_", "", literal=True)
    )

    # Join in the phrog → (product, function) mapping table.
    # ``null_values=["NA"]`` matches pandas' default ``read_csv`` behaviour
    # of interpreting the literal string ``"NA"`` as a missing value, so
    # the subsequent fill_null() below replaces those with the canonical
    # ``"hypothetical protein"`` / ``"unknown function"`` defaults rather
    # than letting the literal ``"NA"`` leak into the output.
    phrog_mapping = pl.read_csv(
        Path(database) / "phold_annots.tsv",
        separator="\t",
        schema_overrides={"phrog": pl.Utf8},
        infer_schema_length=10_000,
        null_values=["NA"],
    )
    foldseek_df = foldseek_df.join(phrog_mapping, on="phrog", how="left").with_columns(
        pl.col("product").fill_null("hypothetical protein"),
        pl.col("function").fill_null("unknown function"),
    )

    # vfdb / card hits get the stricter ``card_vfdb_evalue`` threshold
    # (Enault et al.); all other hits pass through.
    cv_threshold = float(card_vfdb_evalue)
    foldseek_df = foldseek_df.filter(
        (
            (pl.col("phrog").is_in(["vfdb", "card"]))
            & (pl.col("evalue").cast(pl.Float64) < cv_threshold)
        )
        | (~pl.col("phrog").is_in(["vfdb", "card"]))
    )

    # ── pre-compute the row-index so groupby tie-breaking matches pandas ──
    # pandas ``groupby.apply`` + ``idxmin`` returns the FIRST row at the
    # minimum value (stable selection). To get the same behaviour from
    # ``sort + group_by.first`` we need a secondary sort key that mirrors
    # the input row order.
    foldseek_df = foldseek_df.with_row_index("_orig_idx")

    # ── custom_nsmallest equivalent ─────────────────────────────────────────
    # Logic: per query, drop "hypothetical protein" rows unless that's
    # ALL there is. From whatever remains, pick the row with the min
    # evalue (breaking ties on original row order).
    non_hypo = foldseek_df.filter(pl.col("product") != "hypothetical protein")
    non_hypo_queries = non_hypo.select(pl.col("query").unique())
    # Queries where every row was hypothetical — keep all rows for these.
    all_hypo_queries = (
        foldseek_df.select(pl.col("query").unique())
        .join(non_hypo_queries, on="query", how="anti")
    )
    all_hypo_rows = foldseek_df.join(all_hypo_queries, on="query", how="semi")

    candidates = pl.concat([non_hypo, all_hypo_rows])
    topfunction_pl = (
        candidates
        .with_columns(pl.col("evalue").cast(pl.Float64).alias("_evalue_f"))
        .sort(["_evalue_f", "_orig_idx"])
        .group_by("query", maintain_order=False)
        .first()
        .drop(["_evalue_f", "_orig_idx"])
        .sort("query")  # match pandas' default groupby(sort=True) output order
    )

    # Format evalue as scientific 3-decimal-place strings, e.g. "8.203e-50".
    # map_elements is needed because polars has no built-in equivalent of
    # Python's ``"{:.3e}".format``. The dataframe is small (one row per
    # query, typically thousands) so the Python-callback overhead is fine.
    topfunction_pl = topfunction_pl.with_columns(
        pl.col("evalue")
        .cast(pl.Float64)
        .map_elements(lambda v: f"{v:.3e}", return_dtype=pl.Utf8)
    )

    # Drop the row-index column from foldseek_df before downstream use —
    # it was only needed for groupby tie-breaking above.
    foldseek_df = foldseek_df.drop("_orig_idx")

    # ── weighted_function equivalent ────────────────────────────────────────
    # Per query, compute the proportion of total non-unknown-function
    # bitscore claimed by each functional category. Rewritten as a
    # native polars pivot — no per-group Python callback.
    weighted_pl = _build_weighted_bitscore_df(foldseek_df)

    return topfunction_pl, weighted_pl


def _build_weighted_bitscore_df(foldseek_df: pl.DataFrame) -> pl.DataFrame:
    """Polars equivalent of the original ``weighted_function`` groupby.apply.

    Per query, compute:
      - ``<category>_bitscore_proportion`` for each known category — the
        rounded-to-3dp ratio of that function's total bitscore to the
        query's total NON-unknown-function bitscore. Always 0 for
        ``unknown function`` (matches the original logic which skipped
        that key when populating ``weighted_counts_normalised``).
      - ``function_with_highest_bitscore_proportion`` — the
        category-name with the highest proportion (excluding unknown);
        falls back to ``"unknown function"`` when no functional bitscore
        exists.
      - ``top_bitscore_proportion_not_unknown`` — the proportion value
        of that top category (0 when no functional bitscore exists).

    Also reproduces the pandas-specific int-vs-float column-dtype quirk:
    a category column where every value is 0 ends up as int64 (writes as
    ``"0"``); a category column where any query has a nonzero proportion
    ends up as float64 (writes as ``"0.0"`` / ``"1.0"`` / etc.).
    """
    # Sum bitscore per (query, function).
    bitscore_qf = (
        foldseek_df
        .group_by(["query", "function"])
        .agg(pl.col("bitscore").sum().alias("fn_bitscore"))
    )

    # Total non-unknown-function bitscore per query.
    total_functional = (
        bitscore_qf
        .filter(pl.col("function") != "unknown function")
        .group_by("query")
        .agg(pl.col("fn_bitscore").sum().alias("total_functional_bitscore"))
    )

    bp = bitscore_qf.join(total_functional, on="query", how="left").with_columns(
        pl.col("total_functional_bitscore").fill_null(0)
    )

    # Per-(query, function) weighted proportion. Zero when function is
    # "unknown function" or when total_functional_bitscore == 0.
    bp = bp.with_columns(
        pl.when(
            (pl.col("function") == "unknown function")
            | (pl.col("total_functional_bitscore") == 0)
        )
        .then(pl.lit(0.0))
        .otherwise(
            (pl.col("fn_bitscore") / pl.col("total_functional_bitscore")).round(3)
        )
        .alias("weighted_proportion")
    )

    # Top non-unknown function per query (where one exists).
    top_per_query = (
        bp
        .filter(pl.col("function") != "unknown function")
        .filter(pl.col("total_functional_bitscore") > 0)
        .sort(["query", "weighted_proportion"], descending=[False, True])
        .group_by("query", maintain_order=False)
        .first()
        .select([
            "query",
            pl.col("function").alias("function_with_highest_bitscore_proportion"),
            pl.col("weighted_proportion").alias("top_bitscore_proportion_not_unknown"),
        ])
    )

    # Pivot proportions to wide form, one column per category.
    wide = (
        bp.pivot(values="weighted_proportion", index="query", on="function")
        .fill_null(0.0)
    )

    # Anchor on every query that survived the upstream filtering.
    all_queries = foldseek_df.select(pl.col("query").unique()).sort("query")

    result = (
        all_queries
        .join(top_per_query, on="query", how="left")
        .with_columns(
            pl.col("function_with_highest_bitscore_proportion").fill_null(
                "unknown function"
            ),
            pl.col("top_bitscore_proportion_not_unknown").fill_null(0).cast(pl.Float64),
        )
        .join(wide, on="query", how="left")
    )

    # Strip ``.pdb`` / ``.cif`` suffixes from query (matches pandas).
    result = result.with_columns(
        pl.col("query")
        .str.replace_all(".pdb", "", literal=True)
        .str.replace_all(".cif", "", literal=True)
    )

    # ── pandas int/float column-dtype preservation ──────────────────────────
    # pandas concat of per-group 1-row DataFrames produces:
    #   - int64 column when EVERY row was the .get(key, 0) default (== 0 int)
    #   - float64 column when ANY row had a computed proportion (rounded float)
    # to_csv then writes "0" vs "0.0" respectively. The existing snapshots
    # are pinned to that exact formatting. The polars-native rule that
    # gives the same byte-identical output: a category column is int iff
    # the underlying ``function`` value never appeared in the data.
    functions_in_data = set(foldseek_df.select(pl.col("function").unique()).to_series().to_list())

    # Build the final column projection in the canonical order.
    # ``wide`` columns are named after the function values themselves
    # (e.g. "DNA, RNA and nucleotide metabolism"); map each to its long
    # output column name via _BITSCORE_CATEGORIES.
    select_exprs: List[pl.Expr] = [
        pl.col("query"),
        pl.col("function_with_highest_bitscore_proportion"),
        pl.col("top_bitscore_proportion_not_unknown"),
    ]
    cast_to_int: List[str] = []
    for fn_name, col_name in _BITSCORE_CATEGORIES:
        if fn_name in result.columns:
            select_exprs.append(pl.col(fn_name).fill_null(0.0).alias(col_name))
        else:
            # No query had this category — column missing from pivot.
            select_exprs.append(pl.lit(0).alias(col_name))
        if fn_name not in functions_in_data:
            cast_to_int.append(col_name)

    result = result.select(select_exprs)

    if cast_to_int:
        result = result.with_columns(
            *[pl.col(c).cast(pl.Int64) for c in cast_to_int]
        )

    return result


def calculate_topfunctions_results(
    filtered_tophits_df: pl.DataFrame,
    cds_dict: Dict[str, Dict[str, dict]],
    output: Path,
    structures: bool,
    proteins_flag: bool,
    fasta_flag: bool,
) -> Tuple[Dict[str, Dict[str, dict]], pl.DataFrame, Dict[str, Dict[str, str]]]:
    """
    Calculate top function results based on filtered top hits DataFrame
    and update CDS dictionary accordingly.

    The previous ``df.iterrows()`` loop is replaced by
    ``polars.iter_rows(named=True)`` which yields dicts directly. The
    per-CDS lookup against ``filtered_tophits_df`` (previously O(N) inside
    an O(M) loop, i.e. O(N×M)) is replaced by a pre-built dict for O(1)
    access.

    Args:
        filtered_tophits_df (pl.DataFrame): DataFrame containing filtered top hits.
        cds_dict (Dict[str, Dict[str, dict]]): Dictionary containing CDS information.
        output (Path): Output path.
        structures (bool): Indicates whether the input is is in .pdb or .cif format.
        proteins_flag (bool): Indicates whether the input is proteins.
        fasta_flag (bool): Indicates whether the input is in FASTA format.

    Returns:
        (updated_cds_dict, filtered_tophits_df, source_dict)
    """

    # ── build result_dict skeleton from cds_dict (Python only) ────────────
    result_dict: Dict[str, Dict[str, dict]] = {}
    cds_record_dict: Dict[str, str] = {}
    for record_id, cds_entries in cds_dict.items():
        # issue 86 — restore the literal '|' that was mangled upstream
        record_id = record_id.replace("~PIPE~", "|")
        result_dict[record_id] = {}
        for cds_id in cds_entries:
            result_dict[record_id][cds_id] = {}
            cds_record_dict[cds_id] = record_id

    # ── do the structures merge + column reorder in polars ────────────────
    tophits = filtered_tophits_df

    if structures:
        # Add contig_id to each tophit row by joining the cds_id → record_id
        # mapping that we just built.
        cds_record_pl = pl.DataFrame(
            {
                "cds_id": list(cds_record_dict.keys()),
                "contig_id": list(cds_record_dict.values()),
            }
        )
        tophits = tophits.join(cds_record_pl, on="cds_id", how="left")

    # The original pandas code only actually applied the reorder in the
    # non-proteins branch (the proteins-flag column_order was computed but
    # never used to reindex). Preserve that to keep CLI outputs identical.
    if not proteins_flag:
        cols = tophits.columns
        front = [c for c in ("contig_id", "cds_id") if c in cols]
        rest = [c for c in cols if c not in front]
        tophits = tophits.select(front + rest)

    # ── first iter-loop: build result_dict from tophit rows ───────────────
    # ``iter_rows(named=True)`` yields a fresh dict per row (no Series
    # construction overhead like pandas' iterrows).
    for row in tophits.iter_rows(named=True):
        record_id = "proteins" if proteins_flag else row["contig_id"]
        cds_id = row["cds_id"]
        values_dict = {k: row[k] for k in (
            "phrog", "product", "function", "tophit_protein",
            "bitscore", "fident", "evalue",
            "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        )}

        # null record_id → the structure in structure_dir has a foldseek hit
        # but can't be matched up to a contig. polars uses None for nulls.
        if record_id is None:
            if structures:
                logger.warning(
                    f"{cds_id} has a foldseek hit but no record_id was found. "
                    "Please check the way you named the structure file."
                )
            else:  # shouldn't happen for non-structures input
                logger.warning(
                    f"{cds_id} has a foldseek hit but no record_id was found. "
                    "Please check your input."
                )
        else:
            result_dict[record_id][cds_id] = values_dict

    # ── pre-build (contig_id, cds_id) → annotation_source lookup ──────────
    # Used inside the second loop below for O(1) access instead of scanning
    # ``filtered_tophits_df`` once per CDS (the original O(N×M) lookup).
    # Only built when needed; in proteins_flag mode the value is hardcoded
    # to "foldseek" anyway and contig_id may be absent.
    annotation_source_lookup: Dict[Tuple[str, str], str] = {}
    if not proteins_flag and "annotation_source" in tophits.columns:
        annotation_source_lookup = dict(
            zip(
                zip(
                    tophits["contig_id"].to_list(),
                    tophits["cds_id"].to_list(),
                ),
                tophits["annotation_source"].to_list(),
            )
        )

    # Return the (possibly reordered/joined) tophits as polars; the caller
    # will work with it natively.
    filtered_tophits_df = tophits

    # Mutate ``cds_dict`` in place. The previous ``copy.deepcopy(cds_dict)``
    # cloned every BioPython ``SeqFeature`` (including its ``FeatureLocation``
    # and ``qualifiers`` dict) for every CDS — a >1 GB RAM cost on
    # 50k-CDS pangenome inputs and the single biggest avoidable allocation
    # in the compare pipeline. Verified safe because the only caller
    # (``subcommand_compare``) never reads ``cds_dict`` again after this
    # function returns — it works with the returned ``updated_cds_dict``
    # exclusively. Aliasing instead of cloning preserves the rest of the
    # function body unchanged so the diff stays minimal and the original
    # naming (``updated_cds_dict`` reflects intent) is kept.
    updated_cds_dict = cds_dict

    phrog_function_mapping = {
        "unknown function": "unknown function",
        "transcription regulation": "transcription regulation",
        "tail": "tail",
        "other": "other",
        "moron, auxiliary metabolic gene and host takeover": "moron, auxiliary metabolic gene and host takeover",
        "lysis": "lysis",
        "integration and excision": "integration and excision",
        "head and packaging": "head and packaging",
        "DNA, RNA and nucleotide metabolism": "DNA, RNA and nucleotide metabolism",
        "connector": "connector",
    }

    # source dict
    source_dict = {}

    # iterates over the records
    for record_id, record in updated_cds_dict.items():
        source_dict[record_id] = {}
        # iterates over the features
        for cds_id, cds_feature in updated_cds_dict[record_id].items():
            # proteins/FASTA input -> no pharokka input -> fake input to make the updating work
            if proteins_flag is True or fasta_flag is True:
                cds_feature.qualifiers["function"] = ["unknown function"]
                cds_feature.qualifiers["phrog"] = ["No_PHROG"]
                cds_feature.qualifiers["product"] = ["hypothetical protein"]

            # if the result_dict is not empty
            # this is a foldseek hit
            if result_dict[record_id][cds_id] != {}:
                # get the foldseek function
                # function will be None if there is no foldseek hit - shouldn't happen here but error handling
                foldseek_phrog = result_dict[record_id][cds_id].get("phrog", None)

                # add annotation source — proteins-compare is always foldseek;
                # otherwise look up via the pre-built (contig_id, cds_id) dict
                # (O(1) instead of scanning the full df per CDS).
                if proteins_flag:
                    source_dict[record_id][cds_id] = "foldseek"
                else:
                    source_dict[record_id][cds_id] = annotation_source_lookup[
                        (record_id, cds_id)
                    ]

                # same phrog as pharokka - only update the function with new annots in v1 in case users have older pharokka
                if foldseek_phrog == cds_feature.qualifiers["phrog"][0]:
                    updated_cds_dict[record_id][cds_id].qualifiers["product"][0] = (
                        result_dict[record_id][cds_id]["product"]
                    )
                    updated_cds_dict[record_id][cds_id].qualifiers["function"][0] = (
                        result_dict[record_id][cds_id]["function"]
                    )

                # different phrog as pharokka
                else:
                    # where there was no phrog in pharokka
                    if cds_feature.qualifiers["phrog"][0] == "No_PHROG":
                        updated_cds_dict[record_id][cds_id].qualifiers["phrog"][0] = (
                            result_dict[record_id][cds_id]["phrog"]
                        )

                        # Handle missing product qualifier - see https://www.ncbi.nlm.nih.gov/nuccore/OY726582.1/ CAJ1523274.1
                        if "product" not in updated_cds_dict[record_id][cds_id].qualifiers:
                            updated_cds_dict[record_id][cds_id].qualifiers["product"] = [
                                result_dict[record_id][cds_id]["product"]
                            ]
                        else:
                            updated_cds_dict[record_id][cds_id].qualifiers["product"][0] = (
                                result_dict[record_id][cds_id]["product"]
                            )
                        # Handle missing function qualifier
                        if "function" not in updated_cds_dict[record_id][cds_id].qualifiers:
                            updated_cds_dict[record_id][cds_id].qualifiers["function"] = [
                                result_dict[record_id][cds_id]["function"]
                            ]
                        else:
                            updated_cds_dict[record_id][cds_id].qualifiers["function"][0] = (
                                result_dict[record_id][cds_id]["function"]
                            )

                    # pharokka has a phrog (or genbank doesn't exist)
                    else:
                        # if the foldseek result is not unknown function then update
                        if (
                            result_dict[record_id][cds_id]["function"]
                            != "unknown function"
                        ):
                            # if from pharokka input gbk
                            try:
                                # update
                                updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                    0
                                ] = result_dict[record_id][cds_id]["phrog"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "product"
                                ][0] = result_dict[record_id][cds_id]["product"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "function"
                                ][0] = result_dict[record_id][cds_id]["function"]
                            except:  # from Genbank input - won't have phrog or function in the updated_cds_dict. Therefore need to create them
                                updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                    0
                                ] = result_dict[record_id][cds_id]["phrog"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "product"
                                ][0] = result_dict[record_id][cds_id]["product"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "function"
                                ][0] = result_dict[record_id][cds_id]["function"]

                        # if foldseek result has unknown function
                        else:
                            # if the pharokka has unknown function, update with foldseek hit anyway
                            if (
                                cds_feature.qualifiers["function"][0]
                                == "unknown function"
                            ):
                                updated_cds_dict[record_id][cds_id].qualifiers["phrog"][
                                    0
                                ] = result_dict[record_id][cds_id]["phrog"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "product"
                                ][0] = result_dict[record_id][cds_id]["product"]
                                updated_cds_dict[record_id][cds_id].qualifiers[
                                    "function"
                                ][0] = result_dict[record_id][cds_id]["function"]

                            else:
                                # this will be when pharokka does have a known function
                                # keep the pharokka annotation - aka do nothing to the dictionary
                                # but need to update annotation source dict to pharokka
                                source_dict[record_id][cds_id] = "pharokka"

            else:
                # no foldseek hits - empty results dict for the record and cds id
                # so the copy from before will be fine
                # therefore just leave whatever pharokka has (no change to the result_dict)
                # need to add to source dict
                if cds_feature.qualifiers["phrog"][0] == "No_PHROG":
                    source_dict[record_id][cds_id] = "none"
                else:
                    source_dict[record_id][cds_id] = "pharokka"

    return updated_cds_dict, filtered_tophits_df, source_dict


def initialize_function_counts_dict(
    record_id: str, count_dict: Dict[str, int], cds_count: int
) -> Dict[str, int]:
    """
    Initialize function counts dictionary for a given record ID.

    Args:
        record_id (str): ID of the record.
        count_dict (Dict[str, int]): Dictionary containing function counts.
        cds_count (int): Count of CDS.

    Returns:
        Dict[str, int]: Updated function counts dictionary.
    """
    count_dict[record_id]["cds_count"] = cds_count
    count_dict[record_id].update(
        {
            "phrog_count": 0,
            "connector": 0,
            "DNA, RNA and nucleotide metabolism": 0,
            "head and packaging": 0,
            "integration and excision": 0,
            "lysis": 0,
            "moron, auxiliary metabolic gene and host takeover": 0,
            "other": 0,
            "tail": 0,
            "transcription regulation": 0,
            "unknown function": 0,
        }
    )

    return count_dict


#####
# custom db
#####


def get_topcustom_hits(
    result_tsv: Path,
    structures: bool,
    proteins_flag: bool,
) -> pl.DataFrame:
    """Process Foldseek output to extract top custom-DB hits.

    The original
    ``foldseek_df.loc[foldseek_df.groupby("query")["evalue"].idxmin()]``
    pattern is replaced by ``sort.group_by.first()`` with a stable
    secondary sort on row order to reproduce pandas' idxmin tie-break.

    Args:
        result_tsv (Path): Path to the Foldseek custom result TSV file.
        structures (bool): Flag indicating whether structures have been added.
        proteins_flag (bool): Flag indicating whether proteins are used.

    Returns:
        pl.DataFrame: DataFrame containing the top hits extracted from the custom Foldseek output.
    """
    logger.info("Processing Foldseek output")

    base_cols = [
        "query", "target", "bitscore", "fident", "evalue",
        "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
    ]
    col_list = base_cols + (["alntmscore", "lddt"] if structures else [])

    foldseek_df = pl.read_csv(
        result_tsv,
        separator="\t",
        has_header=False,
        new_columns=col_list,
        schema_overrides={"evalue": pl.Utf8},
        infer_schema_length=10_000,
    )

    if foldseek_df.is_empty():
        logger.warning(
            "Foldseek found no custom hits whatsoever - please check your "
            "custom database and input."
        )
        logger.warning("Phold will continue using only the default databases.")

    # issue #86 — restore the literal '|' that was mangled upstream
    foldseek_df = foldseek_df.with_columns(
        pl.col("query").str.replace_all("~PIPE~", "|", literal=True)
    )

    # Derive cds_id from query.
    if not structures and not proteins_flag:
        # prostt5 path: query = "<contig_id>:<cds_id>"
        # NOTE — the original code splits into contig_id+cds_id then
        # IMMEDIATELY drops contig_id. So we just produce cds_id directly.
        foldseek_df = foldseek_df.with_columns(
            pl.col("query")
            .str.splitn(":", 2)
            .struct.rename_fields(["_contig_id", "cds_id"])
            .alias("_q")
        ).unnest("_q").drop("_contig_id")
    else:
        # structures / proteins_flag: query is the cds_id (possibly with
        # ``.pdb``/``.cif`` suffix). Preserve the original pandas bug
        # where the second assignment overwrites the first — only
        # ``.cif`` actually gets stripped from the original query.
        foldseek_df = foldseek_df.with_columns(
            pl.col("query").str.replace_all(".cif", "", literal=True).alias("cds_id")
        )

    # Clean up ``.pdb`` and ``.cif`` suffixes from target.
    foldseek_df = foldseek_df.with_columns(
        pl.col("target")
        .str.replace_all(".pdb", "", literal=True)
        .str.replace_all(".cif", "", literal=True)
    )

    # Pick the lowest-evalue row per query. Add _orig_idx so ties break
    # on input order, matching pandas' idxmin (stable).
    tophit_custom = (
        foldseek_df
        .with_row_index("_orig_idx")
        .with_columns(pl.col("evalue").cast(pl.Float64).alias("_evalue_f"))
        .sort(["_evalue_f", "_orig_idx"])
        .group_by("query", maintain_order=False)
        .first()
        .drop(["_evalue_f", "_orig_idx"])
        .sort("query")  # match pandas groupby(sort=True) output order
    )

    # Drop the query column (downstream uses cds_id instead).
    tophit_custom = tophit_custom.drop("query")

    return tophit_custom


def calculate_qcov_tcov(merged_df: pl.DataFrame) -> pl.DataFrame:
    """Add ``qCov`` and ``tCov`` columns to a foldseek tophits frame.

    ``qCov = round((qEnd - qStart) / qLen, 2)`` and analogously for tCov.
    The columns are inserted immediately after their ``qLen``/``tLen``
    siblings via an explicit reorder.

    Output is byte-identical to the prior pandas implementation
    (verified by
    ``tests/unit/test_topfunction.py::test_calculate_qcov_tcov_snapshot``).
    """
    pl_df = merged_df.with_columns(
        ((pl.col("qEnd") - pl.col("qStart")) / pl.col("qLen")).round(2).alias("qCov"),
        ((pl.col("tEnd") - pl.col("tStart")) / pl.col("tLen")).round(2).alias("tCov"),
    )

    # Reorder: insert qCov directly after qLen, tCov directly after tLen.
    cols = pl_df.columns
    qLen_idx = cols.index("qLen")
    tLen_idx = cols.index("tLen")
    repositioned = {"qCov", "tStart", "tEnd", "tLen", "tCov"}
    before = [c for c in cols[: qLen_idx + 1] if c not in repositioned]
    after  = [c for c in cols[tLen_idx + 1 :] if c not in repositioned]
    new_order = before + ["qCov", "tStart", "tEnd", "tLen", "tCov"] + after

    return pl_df.select(new_order)
