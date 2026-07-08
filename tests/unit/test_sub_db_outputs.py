"""
Function-level tests for phold.io.sub_db_outputs.create_sub_db_outputs.

Five sub-databases (ACR / VFDB / CARD / NetFlax / DefenseFinder) each get
their own TSV written under sub_db_tophits/. Each sub-DB has two paths:
  - hit found -> merge with metadata, write merged TSV
  - no hit    -> touch_file (empty placeholder)

Both paths are gated by snapshot assertions.

NOTE — the relative output path `sub_db_tophits/<name>_cds_predictions.tsv`
is duplicated between production (src/phold/io/sub_db_outputs.py) and this
test file's `_read_subdb_tsv`. The polars migration is a good moment to
extract a `sub_db_output_paths()` helper in the production module and have
both producer and test consume it.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from phold.io.sub_db_outputs import create_sub_db_outputs

SNAPSHOT_SUBDIR = "test_sub_db_outputs"
SUBDBS = ["acr", "vfdb", "card", "netflax", "defensefinder"]


def _read_subdb_tsv(output_dir: Path, name: str) -> str:
    p = output_dir / "sub_db_tophits" / f"{name}_cds_predictions.tsv"
    assert p.exists(), f"expected {p} to exist (even if empty)"
    return p.read_text()


def test_create_sub_db_outputs_returns_true(subdb_input_with_hits_df, tmp_path, phold_db):
    assert create_sub_db_outputs(subdb_input_with_hits_df, phold_db, tmp_path) is True


def test_create_sub_db_outputs_creates_all_five_files(subdb_input_with_hits_df, tmp_path, phold_db):
    """Every sub-DB must produce a file (possibly empty) regardless of input."""
    create_sub_db_outputs(subdb_input_with_hits_df, phold_db, tmp_path)
    for name in SUBDBS:
        path = tmp_path / "sub_db_tophits" / f"{name}_cds_predictions.tsv"
        assert path.exists(), f"missing {path}"


@pytest.mark.parametrize("subdb", SUBDBS)
def test_create_sub_db_outputs_with_hits_snapshot(
    subdb_input_with_hits_df, tmp_path, phold_db, snapshot, subdb
):
    """When the input has one hit per sub-DB, each TSV is non-empty and
    matches its snapshot. This is the byte-identity gate for the polars
    merge/drop chain in each sub-DB block."""
    create_sub_db_outputs(subdb_input_with_hits_df, phold_db, tmp_path)
    snapshot.snapshot_dir = Path(__file__).parent / "snapshots" / SNAPSHOT_SUBDIR
    snapshot.assert_match(_read_subdb_tsv(tmp_path, subdb), f"with_hits_{subdb}.tsv")


@pytest.mark.parametrize("subdb", SUBDBS)
def test_create_sub_db_outputs_no_hits_snapshot(
    subdb_input_no_hits_df, tmp_path, phold_db, snapshot, subdb
):
    """When the input has no hits for any sub-DB (real production case for
    most phages), every sub-DB file must still be created via touch_file.
    Snapshot pins the empty-file convention."""
    create_sub_db_outputs(subdb_input_no_hits_df, phold_db, tmp_path)
    snapshot.snapshot_dir = Path(__file__).parent / "snapshots" / SNAPSHOT_SUBDIR
    snapshot.assert_match(_read_subdb_tsv(tmp_path, subdb), f"no_hits_{subdb}.tsv")
