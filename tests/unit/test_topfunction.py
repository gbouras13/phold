"""
Function-level tests for phold.results.topfunction.

Each test pins one branch of the analytical core via snapshot assertion.
Polars migration will swap implementations; these snapshots are the
byte-identity gate.

Run::

    pytest tests/unit/test_topfunction.py -v
    pytest tests/unit/test_topfunction.py --snapshot-update   # bootstrap / accept changes
"""
from __future__ import annotations

import polars as pl
import pytest

from phold.results.topfunction import (calculate_qcov_tcov, get_topcustom_hits,
                                       get_topfunctions)

SNAPSHOT_SUBDIR = "test_topfunction"


# ─── get_topfunctions ───────────────────────────────────────────────────────


def test_get_topfunctions_basic(fixtures_dir, df_snapshot, phold_db):
    """5 hits across 3 queries; every query has at least one non-hypothetical
    product so the else-branch of custom_nsmallest fires."""
    topfn_df, weighted_df = get_topfunctions(
        fixtures_dir / "foldseek_basic.tsv",
        phold_db,
        database_name="phold",
        structures=False,
        card_vfdb_evalue=1e-10,
        proteins_flag=False,
    )
    df_snapshot(topfn_df,    "basic_topfunction.tsv",       subdir=SNAPSHOT_SUBDIR)
    df_snapshot(weighted_df, "basic_weighted_bitscore.tsv", subdir=SNAPSHOT_SUBDIR)


def test_get_topfunctions_all_hypothetical_group(fixtures_dir, df_snapshot, phold_db):
    """All hits for a query map to unknown-phrog (hypothetical product) —
    custom_nsmallest's `if all(...)` branch must fire and return the
    min-evalue row anyway."""
    topfn_df, _ = get_topfunctions(
        fixtures_dir / "foldseek_all_hypothetical.tsv",
        phold_db,
        database_name="phold",
        structures=False,
        card_vfdb_evalue=1e-10,
        proteins_flag=False,
    )
    df_snapshot(topfn_df, "all_hypothetical_topfunction.tsv", subdir=SNAPSHOT_SUBDIR)


def test_get_topfunctions_vfdb_card_evalue_filter(fixtures_dir, df_snapshot, phold_db):
    """vfdb / card hits with evalue above the stricter card_vfdb_evalue
    threshold must be filtered out before the topfunction selection."""
    topfn_df, weighted_df = get_topfunctions(
        fixtures_dir / "foldseek_vfdb_card_filter.tsv",
        phold_db,
        database_name="phold",
        structures=False,
        card_vfdb_evalue=1e-10,
        proteins_flag=False,
    )
    df_snapshot(topfn_df,    "vfdb_card_filter_topfunction.tsv", subdir=SNAPSHOT_SUBDIR)
    df_snapshot(weighted_df, "vfdb_card_filter_weighted.tsv",    subdir=SNAPSHOT_SUBDIR)


@pytest.fixture(scope="module")
def phrog_prefixes_result(fixtures_dir, phold_db):
    """get_topfunctions on the phrog-prefixes fixture, computed once and
    shared by both prefix-related tests below."""
    topfn_df, _ = get_topfunctions(
        fixtures_dir / "foldseek_phrog_prefixes.tsv",
        phold_db,
        database_name="phold",
        structures=False,
        card_vfdb_evalue=1e-10,
        proteins_flag=False,
    )
    return topfn_df


def test_get_topfunctions_phrog_prefix_stripping(phrog_prefixes_result):
    """envhog_/efam_/dgr_/plain phrog_ prefixes must all reduce to the bare
    numeric phrog id in the output."""
    expected_phrogs = {"453", "4644", "1242"}
    actual = set(phrog_prefixes_result["phrog"].to_list())
    assert actual == expected_phrogs, f"got {actual}"


def test_get_topfunctions_envhog_promotes_to_tophit_protein(phrog_prefixes_result):
    """envhog_ prefix on the target column means the protein name is prefixed
    with `envhog_` in the output's tophit_protein column."""
    # The envhog row had target = 'envhog_phrog_453:1bUP3'.
    # The phrog should be '453' but the tophit_protein should carry the envhog_ marker.
    envhog_proteins = (
        phrog_prefixes_result
        .filter(pl.col("phrog") == "453")
        .get_column("tophit_protein")
        .to_list()
    )
    assert any(p.startswith("envhog_") for p in envhog_proteins), (
        f"expected an envhog_-prefixed tophit_protein in: {envhog_proteins}"
    )


def test_get_topfunctions_pipe_round_trip(fixtures_dir, phold_db):
    """`~PIPE~` markers in the query column must be replaced with `|` early
    so downstream code sees the real contig name."""
    topfn_df, _ = get_topfunctions(
        fixtures_dir / "foldseek_pipe_in_query.tsv",
        phold_db,
        database_name="phold",
        structures=False,
        card_vfdb_evalue=1e-10,
        proteins_flag=False,
    )
    queries = topfn_df["query"].to_list()
    assert all("|" in q for q in queries), f"~PIPE~ not replaced: {queries}"
    assert not any("~PIPE~" in q for q in queries)


# ─── empty Foldseek output (zero hits) ──────────────────────────────────────
# Regression for v1.3.0: Foldseek writes a 0-byte result TSV when nothing
# hits, pl.read_csv raised ``NoDataError: empty CSV`` and phold died. Both
# parsers must return empty frames with the same schema they'd produce for a
# populated TSV, so the rest of compare runs unchanged. Uses mini_phold_db
# (not phold_db) so these run without the database.

# (structures, proteins_flag) — every combination subcommand_compare uses.
_FOLDSEEK_MODES = [(False, False), (False, True), (True, False), (True, True)]
_MODE_IDS = ["prostt5", "proteins", "structures", "proteins_structures"]


def _write_foldseek_tsv(path, structures: bool, with_hit: bool):
    """Write a Foldseek result TSV in the given mode: one hit, or empty."""
    if not with_hit:
        path.write_text("")
        return path
    row = [
        "NC_043029:NC_043029_CDS_0001", "phrog_453:protein221399", "2019",
        "0.934", "8.203e-50", "1", "351", "352", "1", "351", "351",
    ]
    if structures:
        row += ["0.912", "0.854"]  # alntmscore, lddt
    path.write_text("\t".join(row) + "\n")
    return path


@pytest.mark.parametrize("structures,proteins_flag", _FOLDSEEK_MODES, ids=_MODE_IDS)
def test_get_topfunctions_empty_tsv(tmp_path, mini_phold_db, loguru_errors, structures, proteins_flag):
    empty_tsv = _write_foldseek_tsv(tmp_path / "empty.tsv", structures, with_hit=False)
    hit_tsv = _write_foldseek_tsv(tmp_path / "hit.tsv", structures, with_hit=True)

    topfn_df, weighted_df = get_topfunctions(
        empty_tsv, mini_phold_db, "phold", structures, 1e-10, proteins_flag
    )
    ref_topfn_df, ref_weighted_df = get_topfunctions(
        hit_tsv, mini_phold_db, "phold", structures, 1e-10, proteins_flag
    )

    assert topfn_df.is_empty()
    assert weighted_df.is_empty()
    assert topfn_df.schema == ref_topfn_df.schema
    # Category dtypes legitimately differ (Int64 when a function never
    # occurs, see _build_weighted_bitscore_df), so only the columns must match.
    assert weighted_df.columns == ref_weighted_df.columns
    # ERROR exits phold, so zero hits must only warn.
    assert loguru_errors == []


@pytest.mark.parametrize("structures,proteins_flag", _FOLDSEEK_MODES, ids=_MODE_IDS)
def test_get_topcustom_hits_empty_tsv(tmp_path, loguru_errors, structures, proteins_flag):
    empty_tsv = _write_foldseek_tsv(tmp_path / "empty.tsv", structures, with_hit=False)
    hit_tsv = _write_foldseek_tsv(tmp_path / "hit.tsv", structures, with_hit=True)

    custom_df = get_topcustom_hits(empty_tsv, structures, proteins_flag)
    ref_custom_df = get_topcustom_hits(hit_tsv, structures, proteins_flag)

    assert custom_df.is_empty()
    assert custom_df.schema == ref_custom_df.schema
    # compare.py adds coverages to the custom hits next.
    assert calculate_qcov_tcov(custom_df).columns == calculate_qcov_tcov(ref_custom_df).columns
    assert loguru_errors == []


# ─── calculate_qcov_tcov ────────────────────────────────────────────────────


def test_calculate_qcov_tcov_adds_two_columns(qcov_tcov_input_df):
    """The function inserts qCov right after qLen and tCov right after tLen
    via an explicit column reorder. Verify both columns exist and the
    reorder put them in the documented positions:
        ... qLen, qCov, tStart, tEnd, tLen, tCov, ...
    """
    out = calculate_qcov_tcov(qcov_tcov_input_df)
    assert "qCov" in out.columns
    assert "tCov" in out.columns
    cols = list(out.columns)
    assert cols.index("qCov") == cols.index("qLen") + 1, (
        f"qCov not directly after qLen: {cols[cols.index('qLen') - 1: cols.index('qLen') + 3]}"
    )
    assert cols.index("tCov") == cols.index("tLen") + 1, (
        f"tCov not directly after tLen: {cols[cols.index('tLen') - 1: cols.index('tLen') + 3]}"
    )


def test_calculate_qcov_tcov_values(qcov_tcov_input_df):
    """Spot-check the math: qCov = round((qEnd - qStart) / qLen, 2). No +1
    and rounded to two decimal places. tCov analogous."""
    out = calculate_qcov_tcov(qcov_tcov_input_df)
    row_in = qcov_tcov_input_df.row(0, named=True)
    row_out = out.row(0, named=True)
    expected_qcov = round((row_in["qEnd"] - row_in["qStart"]) / row_in["qLen"], 2)
    assert row_out["qCov"] == pytest.approx(expected_qcov, abs=1e-9)
    expected_tcov = round((row_in["tEnd"] - row_in["tStart"]) / row_in["tLen"], 2)
    assert row_out["tCov"] == pytest.approx(expected_tcov, abs=1e-9)


def test_calculate_qcov_tcov_snapshot(qcov_tcov_input_df, df_snapshot):
    """Full snapshot of qcov/tcov output — the byte-identity gate for the
    polars rewrite."""
    out = calculate_qcov_tcov(qcov_tcov_input_df)
    df_snapshot(out, "qcov_tcov_output.tsv", subdir=SNAPSHOT_SUBDIR)
