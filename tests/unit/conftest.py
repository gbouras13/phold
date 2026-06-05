"""
Shared fixtures for the function-level unit tests.

The unit tests are a pre-work scaffold for the planned polars migration:
they pin the current pandas behaviour of every analytical function so any
polars rewrite that drifts shape/dtype/value lights up immediately.

See POLARS_MIGRATION_ANALYSIS.md for the migration plan.
"""
from __future__ import annotations

import io
from pathlib import Path

import polars as pl
import pytest

UNIT_DIR = Path(__file__).parent
FIXTURES = UNIT_DIR / "fixtures"
SNAPSHOTS = UNIT_DIR / "snapshots"
TEST_DATA = UNIT_DIR.parent / "test_data"


# ── path-shaped fixtures ────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def phold_db() -> Path:
    """Path to the test phold database (contains phold_annots.tsv + subdb metadata)."""
    db = TEST_DATA / "phold_db"
    if not db.exists():
        pytest.skip(f"phold test database missing at {db}")
    return db


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    """Directory holding hand-crafted TSV fixtures consumed by tests that
    pass a file path to the function under test (e.g. get_topfunctions)."""
    return FIXTURES


# ── parquet-loaded DataFrame fixtures (session-scoped: read once) ─────────
# Consumers MUST .copy() before passing to mutating functions (e.g.
# calculate_qcov_tcov writes new columns onto the input frame). The
# session scope keeps disk I/O at one read per parquet for the whole run.


@pytest.fixture(scope="session")
def qcov_tcov_input_df() -> pl.DataFrame:
    return pl.read_parquet(FIXTURES / "qcov_tcov_input.parquet")


@pytest.fixture(scope="session")
def subdb_input_with_hits_df() -> pl.DataFrame:
    return pl.read_parquet(FIXTURES / "sub_db_outputs_input_with_all_subdbs.parquet")


@pytest.fixture(scope="session")
def subdb_input_no_hits_df() -> pl.DataFrame:
    return pl.read_parquet(FIXTURES / "sub_db_outputs_input_no_subdbs.parquet")


# ── snapshot helper: canonical TSV + per-module snapshot dir + assertion ──


@pytest.fixture
def df_snapshot(snapshot):
    """Snapshot-compare a polars DataFrame as canonical TSV.

    The canonical form (no row index, tab-separated) is the single
    project-wide convention for DataFrame snapshots so all test files
    round-trip the same way. Snapshots live under
    ``tests/unit/snapshots/<subdir>/<name>``.

    Usage::

        def test_thing(df_snapshot):
            out_df = some_function(...)
            df_snapshot(out_df, "expected.tsv", subdir="test_thing_module")
    """
    def _snapshot(df: pl.DataFrame, name: str, subdir: str) -> None:
        snapshot.snapshot_dir = SNAPSHOTS / subdir
        buf = io.BytesIO()
        df.write_csv(buf, separator="\t")
        tsv = buf.getvalue().decode("utf-8")
        snapshot.assert_match(tsv, name)

    return _snapshot
