"""
Function-level tests for phold.io.tbl — the NCBI feature table writer (issue #137).

The .tbl is a submission format: table2asn parses it positionally, so the
tests below pin the exact line shape (tab counts, coordinate order, which
qualifiers are emitted) rather than just "a file was written".

Run::

    pytest tests/unit/test_tbl.py -v
"""
from __future__ import annotations

import polars as pl
import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature

from phold.io.tbl import write_tbl


def _per_cds_df(rows):
    return pl.DataFrame(
        rows,
        schema={
            "contig_id": pl.Utf8,
            "cds_id": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "strand": pl.Int64,
            "product": pl.Utf8,
            "function": pl.Utf8,
            "annotation_method": pl.Utf8,
            "transl_table": pl.Utf8,
        },
    )


def _cds_row(**kw):
    row = {
        "contig_id": "contig1",
        "cds_id": "contig1_CDS_0001",
        "start": 100,
        "end": 400,
        "strand": 1,
        "product": "terminase large subunit",
        "function": "head and packaging",
        "annotation_method": "foldseek",
        "transl_table": "11",
    }
    row.update(kw)
    return row


def _write(tmp_path, cds_rows, non_cds=None, contig_ids=("contig1",)):
    out = write_tbl(
        _per_cds_df(cds_rows),
        non_cds or {},
        list(contig_ids),
        "phold",
        tmp_path,
    )
    return out.read_text()


# ─── basic shape ────────────────────────────────────────────────────────────


def test_writes_prefixed_filename(tmp_path):
    write_tbl(_per_cds_df([_cds_row()]), {}, ["contig1"], "myprefix", tmp_path)
    assert (tmp_path / "myprefix.tbl").exists()


def test_cds_block_format(tmp_path):
    """One CDS -> a >Feature header plus a coordinate line and 4 qualifiers."""
    text = _write(tmp_path, [_cds_row()])

    assert text == (
        ">Feature contig1\n"
        "100\t400\tCDS\n"
        "\t\t\tproduct\tterminase large subunit\n"
        "\t\t\tfunction\thead and packaging\n"
        "\t\t\tinference\tfoldseek\n"
        "\t\t\ttransl_table\t11\n"
    )


def test_qualifier_lines_have_exactly_three_leading_tabs(tmp_path):
    """table2asn is positional — the qualifier indent must be 3 tabs, not 2 or 4."""
    text = _write(tmp_path, [_cds_row()])
    for line in text.splitlines():
        if line.startswith("\t"):
            assert line.startswith("\t\t\t") and not line.startswith("\t\t\t\t")


# ─── strand encoding ────────────────────────────────────────────────────────


def test_minus_strand_swaps_coordinates(tmp_path):
    """Strand is encoded by coordinate order — minus strand is written high-first."""
    text = _write(tmp_path, [_cds_row(strand=-1, start=100, end=400)])
    assert "400\t100\tCDS\n" in text


def test_plus_strand_keeps_coordinate_order(tmp_path):
    text = _write(tmp_path, [_cds_row(strand=1, start=100, end=400)])
    assert "100\t400\tCDS\n" in text


# ─── null / empty handling ──────────────────────────────────────────────────


def test_null_product_is_skipped_not_written_empty(tmp_path):
    """A bare 'product' line with no value is a table2asn parse error."""
    text = _write(tmp_path, [_cds_row(product=None)])
    assert "\t\t\tproduct\n" not in text
    assert "\t\t\tproduct\t\n" not in text
    assert "\t\t\tfunction\thead and packaging\n" in text


def test_empty_dataframe_still_writes_contig_headers(tmp_path):
    """Contigs with no features still need a >Feature header for table2asn."""
    out = write_tbl(
        pl.DataFrame(schema={"contig_id": pl.Utf8}),
        {},
        ["contig1", "contig2"],
        "phold",
        tmp_path,
    )
    assert out.read_text() == ">Feature contig1\n>Feature contig2\n"


# ─── multi-contig ───────────────────────────────────────────────────────────


def test_features_are_grouped_under_their_own_contig(tmp_path):
    text = _write(
        tmp_path,
        [
            _cds_row(contig_id="contig1", start=1, end=99),
            _cds_row(contig_id="contig2", start=5, end=50),
        ],
        contig_ids=("contig1", "contig2"),
    )
    c1, c2 = text.split(">Feature contig2\n")
    assert "1\t99\tCDS\n" in c1
    assert "5\t50\tCDS\n" in c2
    assert "1\t99\tCDS\n" not in c2


def test_contig_order_follows_supplied_order(tmp_path):
    text = _write(
        tmp_path,
        [_cds_row(contig_id="contig2"), _cds_row(contig_id="contig1")],
        contig_ids=("contig2", "contig1"),
    )
    assert text.index(">Feature contig2") < text.index(">Feature contig1")


# ─── non-CDS pass-through (tRNA / CRISPR / tmRNA) ───────────────────────────


def _trna_feature():
    """A tRNA shaped like the ones in tests/test_data/SAOMS1.gbk."""
    return SeqFeature(
        FeatureLocation(115193, 115265, strand=-1),
        type="tRNA",
        qualifiers={
            "ID": ["DQEPQRKE_tRNA_0001"],
            "transl_table": ["11"],
            "product": ["tRNA-Met(CAT)"],
            "isotype": ["Met"],
            "anticodon": ["CAT"],
            "locus_tag": ["DQEPQRKE_tRNA_0001"],
            "source": ["tRNAscan-SE_2.0.12"],
            "score": ["61.2"],
        },
    )


def test_trna_is_emitted_with_biopython_coordinate_conversion(tmp_path):
    """BioPython locations are 0-based half-open; the .tbl is 1-based inclusive.

    FeatureLocation(115193, 115265) is bases 115194..115265, and on the minus
    strand that is written high-coordinate first.
    """
    text = _write(tmp_path, [], {"contig1": {"t1": _trna_feature()}})
    assert "115265\t115194\ttRNA\n" in text


def test_trna_qualifiers_filtered_and_ordered(tmp_path):
    text = _write(tmp_path, [], {"contig1": {"t1": _trna_feature()}})

    assert "\t\t\tproduct\ttRNA-Met(CAT)\n" in text
    assert "\t\t\tanticodon\tCAT\n" in text
    # product must precede anticodon per _NON_CDS_QUALIFIER_ORDER
    assert text.index("product\ttRNA") < text.index("anticodon\tCAT")


@pytest.mark.parametrize(
    "blocked", ["ID", "locus_tag", "score", "source", "isotype", "transl_table"]
)
def test_internal_qualifiers_never_reach_the_submission_file(tmp_path, blocked):
    """Pipeline bookkeeping must not leak into a GenBank submission."""
    text = _write(tmp_path, [], {"contig1": {"t1": _trna_feature()}})
    assert f"\t\t\t{blocked}\t" not in text


def test_biopython_bare_string_qualifiers_are_handled(tmp_path):
    """Phold-built features store qualifiers as bare strings, not lists."""
    feature = SeqFeature(
        FeatureLocation(0, 60, strand=1),
        type="tRNA",
        qualifiers={"product": "tRNA-Phe(GAA)"},
    )
    text = _write(tmp_path, [], {"contig1": {"t": feature}})
    assert "1\t60\ttRNA\n" in text
    assert "\t\t\tproduct\ttRNA-Phe(GAA)\n" in text


def test_pharokka_trna_qualifier_becomes_product(tmp_path):
    """Pharokka's /trna= is not a valid NCBI qualifier; product is, and is
    mandatory for tRNA. write_genbank renames it in place, but write_tbl must
    not depend on that having run first."""
    feature = SeqFeature(
        FeatureLocation(0, 72, strand=1),
        type="tRNA",
        qualifiers={"trna": ["tRNA-Met(CAT)"], "anticodon": ["CAT"]},
    )
    text = _write(tmp_path, [], {"contig1": {"t": feature}})

    assert "\t\t\tproduct\ttRNA-Met(CAT)\n" in text
    assert "\t\t\ttrna\t" not in text


def test_existing_product_wins_over_trna_alias(tmp_path):
    feature = SeqFeature(
        FeatureLocation(0, 72, strand=1),
        type="tRNA",
        qualifiers={"trna": ["stale"], "product": ["tRNA-Phe(GAA)"]},
    )
    text = _write(tmp_path, [], {"contig1": {"t": feature}})

    assert "\t\t\tproduct\ttRNA-Phe(GAA)\n" in text
    assert "stale" not in text


def test_write_tbl_does_not_mutate_input_features(tmp_path):
    """The caller's SeqFeature dict is shared with the GenBank writer."""
    feature = SeqFeature(
        FeatureLocation(0, 72, strand=1),
        type="tRNA",
        qualifiers={"trna": ["tRNA-Met(CAT)"]},
    )
    _write(tmp_path, [], {"contig1": {"t": feature}})

    assert feature.qualifiers == {"trna": ["tRNA-Met(CAT)"]}


def test_crispr_repeat_region_qualifiers(tmp_path):
    feature = SeqFeature(
        FeatureLocation(999, 1050, strand=1),
        type="repeat_region",
        qualifiers={
            "rpt_family": ["CRISPR"],
            "rpt_type": ["direct"],
            "rpt_unit_seq": ["GTTTC"],
            "ID": ["crispr_1"],
        },
    )
    text = _write(tmp_path, [], {"contig1": {"c": feature}})
    assert "1000\t1050\trepeat_region\n" in text
    assert "\t\t\trpt_family\tCRISPR\n" in text
    assert "\t\t\trpt_unit_seq\tGTTTC\n" in text
    assert "\t\t\tID\t" not in text


def test_non_cds_features_sorted_by_start(tmp_path):
    late = SeqFeature(FeatureLocation(500, 560, strand=1), type="tRNA", qualifiers={})
    early = SeqFeature(FeatureLocation(10, 70, strand=1), type="tRNA", qualifiers={})
    text = _write(tmp_path, [], {"contig1": {"late": late, "early": early}})
    assert text.index("11\t70\ttRNA") < text.index("501\t560\ttRNA")


def test_cds_and_non_cds_both_present(tmp_path):
    text = _write(tmp_path, [_cds_row()], {"contig1": {"t1": _trna_feature()}})
    assert "100\t400\tCDS\n" in text
    assert "115265\t115194\ttRNA\n" in text
    assert text.count(">Feature contig1\n") == 1
