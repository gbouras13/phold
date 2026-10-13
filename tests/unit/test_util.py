"""
Function-level tests for the file helpers in phold.utils.util.

``read_headerless_csv`` exists because polars can't read an empty header-less
file with ``new_columns`` (``NoDataError``, or ``ShapeError`` with
``raise_if_empty=False``), which crashed phold v1.3.0 whenever Foldseek found
no hits.
"""
from __future__ import annotations

import polars as pl
import pytest

from phold.utils.util import is_blank_file, read_headerless_csv

SCHEMA = {"cds_id": pl.Utf8, "prostt5_confidence": pl.Float64}


@pytest.mark.parametrize("content", ["", "\n", "  \n\t\n"], ids=["empty", "newline", "whitespace"])
def test_read_headerless_csv_blank_file_returns_typed_empty_frame(tmp_path, content):
    path = tmp_path / "blank.csv"
    path.write_text(content)

    assert is_blank_file(path)
    df = read_headerless_csv(path, SCHEMA)

    assert df.is_empty()
    assert dict(df.schema) == SCHEMA


def test_read_headerless_csv_populated_file(tmp_path):
    path = tmp_path / "probs.csv"
    path.write_text("CDS_0001,55.5\nCDS_0002,71.25\n")

    assert not is_blank_file(path)
    df = read_headerless_csv(path, SCHEMA)

    assert dict(df.schema) == SCHEMA
    assert df.rows() == [("CDS_0001", 55.5), ("CDS_0002", 71.25)]


def test_read_headerless_csv_passes_read_csv_kwargs(tmp_path):
    """Extra kwargs (here ``schema_overrides``) reach ``pl.read_csv``."""
    path = tmp_path / "hits.tsv"
    path.write_text("q1\t1e-30\n")

    df = read_headerless_csv(
        path,
        {"query": pl.Utf8, "evalue": pl.Utf8},
        separator="\t",
        schema_overrides={"evalue": pl.Utf8},
    )

    assert df.rows() == [("q1", "1e-30")]
