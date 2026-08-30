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
            # write_genbank converts strand to "+"/"-" strings before the .tbl
            # writer ever sees per_cds_df, so the fixture matches that.
            "strand": pl.Utf8,
            "product": pl.Utf8,
            "function": pl.Utf8,
            "annotation_method": pl.Utf8,
            "transl_table": pl.Utf8,
            "partial": pl.Utf8,
            "phrog": pl.Utf8,
        },
    )


def _cds_row(**kw):
    row = {
        "contig_id": "contig1",
        "cds_id": "contig1_CDS_0001",
        "start": 100,
        "end": 400,
        "strand": "+",
        "product": "terminase large subunit",
        "function": "head and packaging",
        "annotation_method": "foldseek",
        "transl_table": "11",
        "partial": "00",
        "phrog": "1215",
    }
    row.update(kw)
    return row


def _write(
    tmp_path,
    cds_rows,
    non_cds=None,
    contig_ids=("contig1",),
    contig_lengths=None,
):
    out = write_tbl(
        _per_cds_df(cds_rows),
        non_cds or {},
        list(contig_ids),
        "phold",
        tmp_path,
        contig_lengths=contig_lengths,
    )
    return out.read_text()


# ─── basic shape ────────────────────────────────────────────────────────────


def test_writes_prefixed_filename(tmp_path):
    write_tbl(_per_cds_df([_cds_row()]), {}, ["contig1"], "myprefix", tmp_path)
    assert (tmp_path / "myprefix.tbl").exists()


def test_cds_block_format(tmp_path):
    """One CDS -> a >Feature header plus a coordinate line and 4 qualifiers."""
    text = _write(tmp_path, [_cds_row(phrog="1215")])

    assert text == (
        ">Feature contig1\n"
        "100\t400\tCDS\n"
        "\t\t\tproduct\tterminase large subunit\n"
        "\t\t\tfunction\thead and packaging\n"
        "\t\t\tinference\tprotein motif:PHROG:1215\n"
        "\t\t\ttransl_table\t11\n"
    )


def test_qualifier_lines_have_exactly_three_leading_tabs(tmp_path):
    """table2asn is positional — the qualifier indent must be 3 tabs, not 2 or 4."""
    text = _write(tmp_path, [_cds_row()])
    for line in text.splitlines():
        if line.startswith("\t"):
            assert line.startswith("\t\t\t") and not line.startswith("\t\t\t\t")


# ─── strand encoding ────────────────────────────────────────────────────────


def test_minus_strand_is_written_high_coordinate_first(tmp_path):
    """Strand is encoded by coordinate order — minus strand is written high-first.

    per_cds_df is ALREADY in transcription order: write_genbank assigns
    ``start = location.end`` for a minus-strand CDS, so start > end in the
    dataframe and the writer must pass the pair straight through. Re-orienting
    here would write every reverse-strand CDS backwards.
    """
    text = _write(tmp_path, [_cds_row(strand="-", start=400, end=100)])
    assert "400\t100\tCDS\n" in text


def test_plus_strand_keeps_coordinate_order(tmp_path):
    text = _write(tmp_path, [_cds_row(strand="+", start=100, end=400)])
    assert "100\t400\tCDS\n" in text


# ─── partial CDS (incomplete ends) ──────────────────────────────────────────
#
# Prodigal's ``partial`` flag is two digits in *genomic* orientation, left edge
# then right edge, on both strands. Verified against pyrodigal's own GFF
# writer: reverse-complementing a sequence turns a left-edge "10" gene into a
# right-edge "01" gene.


def test_complete_cds_has_no_partial_markers(tmp_path):
    text = _write(tmp_path, [_cds_row(partial="00")])
    assert "100\t400\tCDS\n" in text
    assert "<" not in text and ">Feature" in text


def test_missing_partial_column_is_treated_as_complete(tmp_path):
    """Older Pharokka / NCBI / Bakta input carries no partial qualifier."""
    text = _write(tmp_path, [_cds_row(partial=None)])
    assert "100\t400\tCDS\n" in text
    assert "\t\t\tcodon_start\t" not in text


def test_plus_strand_five_prime_partial(tmp_path):
    """+ strand, partial at the genomic left = incomplete 5' end -> '<'."""
    text = _write(tmp_path, [_cds_row(strand="+", start=1, end=400, partial="10")])
    assert "<1\t400\tCDS\n" in text


def test_plus_strand_three_prime_partial(tmp_path):
    """+ strand, partial at the genomic right = incomplete 3' end -> '>'."""
    text = _write(tmp_path, [_cds_row(strand="+", start=100, end=400, partial="01")])
    assert "100\t>400\tCDS\n" in text


def test_minus_strand_five_prime_partial(tmp_path):
    """- strand: the 5' end is the genomic RIGHT edge, so "01" marks it.

    per_cds_df stores minus-strand CDS as start=high, end=low.
    """
    text = _write(
        tmp_path,
        [_cds_row(strand="-", start=400, end=100, partial="01")],
        contig_lengths={"contig1": 400},
    )
    # the high coordinate is the 5' end and carries the '<'
    assert "<400\t100\tCDS\n" in text


def test_minus_strand_three_prime_partial(tmp_path):
    """- strand: the 3' end is the genomic LEFT edge, so "10" marks it."""
    text = _write(tmp_path, [_cds_row(strand="-", start=400, end=1, partial="10")])
    assert "400\t>1\tCDS\n" in text


def test_both_ends_partial_gets_both_markers(tmp_path):
    """Pharokka's if/elif chain misses "11" entirely; both ends must be marked."""
    text = _write(tmp_path, [_cds_row(strand="+", start=1, end=400, partial="11")])
    assert "<1\t>400\tCDS\n" in text


def test_codon_start_when_plus_partial_not_flush_with_contig_start(tmp_path):
    """A 5'-partial CDS starting at base 2 is extended to <1 with codon_start 2."""
    text = _write(tmp_path, [_cds_row(strand="+", start=2, end=400, partial="10")])
    assert "<1\t400\tCDS\n" in text
    assert "\t\t\tcodon_start\t2\n" in text


def test_no_codon_start_when_plus_partial_flush_at_base_one(tmp_path):
    text = _write(tmp_path, [_cds_row(strand="+", start=1, end=400, partial="10")])
    assert "\t\t\tcodon_start\t" not in text


def test_codon_start_when_minus_partial_not_flush_with_contig_end(tmp_path):
    """Minus-strand mirror: 5' end at 398 on a 400bp contig -> codon_start 3."""
    text = _write(
        tmp_path,
        [_cds_row(strand="-", start=398, end=100, partial="01")],
        contig_lengths={"contig1": 400},
    )
    assert "<400\t100\tCDS\n" in text
    assert "\t\t\tcodon_start\t3\n" in text


def test_minus_partial_without_contig_length_still_marks_but_omits_codon_start(
    tmp_path,
):
    """contig_lengths is optional — degrade to the '<' marker alone."""
    text = _write(tmp_path, [_cds_row(strand="-", start=398, end=100, partial="01")])
    assert "<398\t100\tCDS\n" in text
    assert "\t\t\tcodon_start\t" not in text


def test_minus_and_plus_partials_are_mirror_images(tmp_path):
    """The same gene reverse-complemented must produce mirrored output.

    A + strand 5'-partial CDS at 2..556 on a 7148bp contig gives "<1 556" with
    codon_start 2. Its revcomp is a - strand 5'-partial CDS at 7147..6593,
    which must give "<7148 6593" with the same codon_start.
    """
    plus = _write(
        tmp_path,
        [_cds_row(strand="+", start=2, end=556, partial="10")],
        contig_lengths={"contig1": 7148},
    )
    minus = _write(
        tmp_path,
        [_cds_row(strand="-", start=7147, end=6593, partial="01")],
        contig_lengths={"contig1": 7148},
    )

    assert "<1\t556\tCDS\n" in plus
    assert "<7148\t6593\tCDS\n" in minus
    assert "\t\t\tcodon_start\t2\n" in plus
    assert "\t\t\tcodon_start\t2\n" in minus


def test_codon_start_over_three_still_writes(tmp_path):
    """codon_start is a frame offset; >3 is invalid and is warned about, but
    must not abort at the end of a long run."""
    text = _write(tmp_path, [_cds_row(strand="+", start=9, end=400, partial="10")])
    assert "<1\t400\tCDS\n" in text
    assert "\t\t\tcodon_start\t9\n" in text


def test_codon_start_written_after_transl_table(tmp_path):
    text = _write(tmp_path, [_cds_row(strand="+", start=2, end=400, partial="10")])
    assert text.index("transl_table") < text.index("codon_start")


@pytest.mark.parametrize("bogus", ["", "1", "abc", "123", "0"])
def test_malformed_partial_is_ignored(tmp_path, bogus):
    """A junk partial value must not produce junk coordinates."""
    text = _write(tmp_path, [_cds_row(partial=bogus)])
    assert "100\t400\tCDS\n" in text


def test_three_prime_partial_extends_to_contig_end(tmp_path):
    """NCBI wants a 3'-partial feature to reach the sequence edge.

    Gene callers stop at the last whole codon, leaving trailing bases, and
    table2asn then warns PartialProblemNotSpliceConsensus3Prime.
    """
    text = _write(
        tmp_path,
        [_cds_row(strand="+", start=6965, end=7147, partial="01")],
        contig_lengths={"contig1": 7148},
    )
    assert "6965\t>7148\tCDS\n" in text


def test_three_prime_partial_extends_to_base_one_on_minus_strand(tmp_path):
    text = _write(
        tmp_path,
        [_cds_row(strand="-", start=184, end=2, partial="10")],
        contig_lengths={"contig1": 7148},
    )
    assert "184\t>1\tCDS\n" in text


def test_three_prime_partial_without_contig_length_is_left_alone(tmp_path):
    text = _write(tmp_path, [_cds_row(strand="+", start=100, end=400, partial="01")])
    assert "100\t>400\tCDS\n" in text


# ─── /inference (INSDC controlled vocabulary) ───────────────────────────────
#
# table2asn rejects any /inference without a category prefix from the INSDC
# list, so Phold's raw annotation_method values can never validate on their own.


def test_inference_carries_insdc_category_and_phrog(tmp_path):
    text = _write(tmp_path, [_cds_row(annotation_method="foldseek", phrog="1215")])
    assert "\t\t\tinference\tprotein motif:PHROG:1215\n" in text


def test_inference_omitted_without_a_phrog_to_cite(tmp_path):
    """No accession to point at, and table2asn flags a bare category too."""
    text = _write(tmp_path, [_cds_row(annotation_method="foldseek", phrog="No_PHROG")])
    assert "\t\t\tinference\t" not in text
    assert "100\t400\tCDS\n" in text


def test_inference_omitted_when_nothing_was_annotated(tmp_path):
    """annotation_method 'none' means no evidence can honestly be cited."""
    text = _write(tmp_path, [_cds_row(annotation_method="none", phrog="No_PHROG")])
    assert "\t\t\tinference\t" not in text
    assert "100\t400\tCDS\n" in text


@pytest.mark.parametrize("method", ["foldseek", "pharokka", "card", "vfdb"])
def test_no_raw_method_name_ever_reaches_the_inference_qualifier(tmp_path, method):
    text = _write(tmp_path, [_cds_row(annotation_method=method)])
    assert f"\t\t\tinference\t{method}\n" not in text


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
    # a bare 3-letter anticodon is not a valid NCBI location and is dropped
    assert "\t\t\tanticodon\t" not in text


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


def test_valid_anticodon_location_is_preserved(tmp_path):
    """Pharokka >= v1.8.2 writes the proper NCBI location form, which is kept."""
    feature = SeqFeature(
        FeatureLocation(115193, 115265, strand=-1),
        type="tRNA",
        qualifiers={
            "product": ["tRNA-Met(CAT)"],
            "anticodon": ["(pos:115194..115196,aa:Met,seq:cat)"],
        },
    )
    text = _write(tmp_path, [], {"contig1": {"t": feature}})
    assert "\t\t\tanticodon\t(pos:115194..115196,aa:Met,seq:cat)\n" in text


@pytest.mark.parametrize("bare", ["CAT", "GAA", "gtc"])
def test_bare_anticodon_code_is_dropped(tmp_path, bare):
    """Pharokka < v1.8.2 wrote a bare 3-letter code, which table2asn rejects."""
    feature = SeqFeature(
        FeatureLocation(0, 72, strand=1),
        type="tRNA",
        qualifiers={"product": ["tRNA-Met(CAT)"], "anticodon": [bare]},
    )
    text = _write(tmp_path, [], {"contig1": {"t": feature}})
    assert "\t\t\tanticodon\t" not in text
    # the tRNA itself is still emitted, with its mandatory product
    assert "\t\t\tproduct\ttRNA-Met(CAT)\n" in text


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
