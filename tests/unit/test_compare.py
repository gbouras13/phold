"""
Function-level tests for phold.subcommands.compare.

Focus: the row-wise `assign_annotation_confidence` function. This is the
single place in compare.py that does per-row classification with 8 nested
boolean conditions across two `structures` paths — the highest-risk
single function in the polars migration.

Polars rewrite will turn this into a chain of `pl.when().then()` — these
parametrised tests pin every branch so any precedence error in the
rewrite breaks a named test.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from phold.io.handle_genbank import get_genbank
from phold.subcommands.compare import (
    assign_annotation_confidence,
    assign_annotation_confidence_expr,
    subcommand_compare,
)


def _row(**kwargs):
    """Build a synthetic row dict with defaults for unused fields."""
    base = dict(
        annotation_method="foldseek",
        qCov=0.5,
        tCov=0.5,
        fident=0.2,
        evalue="1e-2",
        prostt5_confidence=20,
    )
    base.update(kwargs)
    return base


# ─── early-return paths (no structures/no structures distinction) ───────────


@pytest.mark.parametrize("structures", [True, False])
def test_annotation_method_none_returns_none(structures):
    assert assign_annotation_confidence(_row(annotation_method="none"), structures=structures) == "none"


@pytest.mark.parametrize("structures", [True, False])
def test_annotation_method_pharokka_returns_pharokka(structures):
    assert assign_annotation_confidence(_row(annotation_method="pharokka"), structures=structures) == "pharokka"


# ─── structures path ────────────────────────────────────────────────────────


def test_structures_high_via_fident():
    """Both coverages > 0.8 AND fident > 0.3 → high."""
    r = _row(qCov=0.9, tCov=0.9, fident=0.5)
    assert assign_annotation_confidence(r, structures=True) == "high"


def test_structures_high_via_evalue():
    """Both coverages > 0.8 AND evalue < 1e-10 (low fident) → high."""
    r = _row(qCov=0.9, tCov=0.9, fident=0.0, evalue="1e-20")
    assert assign_annotation_confidence(r, structures=True) == "high"


def test_structures_medium_via_one_coverage_and_fident():
    """One coverage > 0.8, fident > 0.3 → medium (high needs BOTH coverages)."""
    r = _row(qCov=0.9, tCov=0.5, fident=0.5)
    assert assign_annotation_confidence(r, structures=True) == "medium"


def test_structures_medium_via_one_coverage_and_evalue():
    """One coverage > 0.8, evalue < 1e-5 (but ≥ 1e-10) → medium."""
    r = _row(qCov=0.9, tCov=0.5, fident=0.0, evalue="1e-7")
    assert assign_annotation_confidence(r, structures=True) == "medium"


def test_structures_low_default():
    """Low coverage AND low fident AND weak evalue → low."""
    r = _row(qCov=0.1, tCov=0.1, fident=0.1, evalue="1e-2")
    assert assign_annotation_confidence(r, structures=True) == "low"


# ─── no-structures path (uses prostt5_confidence) ──────────────────────────


def test_no_structures_high_via_prostt5():
    """Both coverages > 0.8 AND prostt5_confidence > 60 → high."""
    r = _row(qCov=0.9, tCov=0.9, fident=0.0, evalue="1e-3", prostt5_confidence=70)
    assert assign_annotation_confidence(r, structures=False) == "high"


def test_no_structures_medium_via_prostt5_in_range_and_evalue():
    """One coverage > 0.8, prostt5 in 45-60, evalue < 1e-5 → medium."""
    r = _row(qCov=0.9, tCov=0.5, fident=0.0, evalue="1e-7", prostt5_confidence=50)
    assert assign_annotation_confidence(r, structures=False) == "medium"


def test_no_structures_low_when_evalue_too_weak():
    """Even with high prostt5, if evalue ≥ 1e-5 then medium tier requires fident — falls to low."""
    r = _row(qCov=0.9, tCov=0.5, fident=0.0, evalue="1e-3", prostt5_confidence=50)
    assert assign_annotation_confidence(r, structures=False) == "low"


def test_no_structures_low_default():
    r = _row(qCov=0.1, tCov=0.1, fident=0.1, evalue="1e-2", prostt5_confidence=10)
    assert assign_annotation_confidence(r, structures=False) == "low"


# ─── parametrised cartesian product for coverage of edge cases ──────────────


@pytest.mark.parametrize(
    "qcov,tcov,fident,evalue,prostt5,structures,expected",
    [
        # ─── structures=True ───────────────────────────────────────────────
        # equal-boundary cases (>, not >=)
        (0.8, 0.9, 0.5, "1e-20", 0, True,  "medium"),   # qCov == 0.8 not > 0.8
        (0.9, 0.8, 0.5, "1e-20", 0, True,  "medium"),   # tCov == 0.8 not > 0.8
        (0.81, 0.81, 0.31, "1", 0, True,  "high"),       # bare-minimum high via fident
        (0.81, 0.81, 0.3,  "1e-11", 0, True, "high"),    # fident exactly at boundary (>0.3 fails) → falls to evalue<1e-10
        (0.81, 0.81, 0.3,  "1e-10", 0, True, "medium"),  # evalue exactly at 1e-10 (>=, not <) → high fails, fident<=0.3 too, but evalue<1e-5 in medium

        # ─── structures=False ──────────────────────────────────────────────
        (0.9, 0.9, 0.0, "1",     61, False, "high"),     # prostt5 > 60
        (0.9, 0.9, 0.0, "1",     60, False, "low"),      # prostt5 == 60 not > 60; no high; no medium (single-cov condition false; evalue weak)
        (0.9, 0.5, 0.0, "1e-7",  45, False, "medium"),   # one cov, prostt5 in range, evalue ok
        (0.9, 0.5, 0.0, "1e-7",  44, False, "low"),      # one cov, prostt5 just below range
        (0.9, 0.5, 0.0, "1e-7",  61, False, "low"),      # one cov, prostt5 above range → medium needs 45..60
    ],
)
def test_confidence_parametrised_boundaries(qcov, tcov, fident, evalue, prostt5, structures, expected):
    """Boundary cases — pin behaviour at the exact > / >= edges."""
    r = _row(qCov=qcov, tCov=tcov, fident=fident, evalue=evalue, prostt5_confidence=prostt5)
    assert assign_annotation_confidence(r, structures=structures) == expected, \
        f"row={r} structures={structures}"


# ─── scalar ⇔ vectorised polars agreement ───────────────────────────────────
# The scalar function ``assign_annotation_confidence`` is exercised by all the
# tests above. ``assign_annotation_confidence_expr`` is the polars vectorised
# rewrite used in production (subcommand_compare). The test below stacks every
# scalar test case into a single DataFrame and confirms the vectorised
# expression returns exactly the same label per row — proving the row-wise
# pandas .apply call we just replaced is semantically equivalent.

# All the rows the scalar tests use, plus their expected labels.
_AGREEMENT_ROWS: list[tuple[dict, bool, str]] = [
    # early-return paths
    (_row(annotation_method="none"),     True,  "none"),
    (_row(annotation_method="none"),     False, "none"),
    (_row(annotation_method="pharokka"), True,  "pharokka"),
    (_row(annotation_method="pharokka"), False, "pharokka"),
    # structures path
    (_row(qCov=0.9, tCov=0.9, fident=0.5),                                    True,  "high"),
    (_row(qCov=0.9, tCov=0.9, fident=0.0, evalue="1e-20"),                    True,  "high"),
    (_row(qCov=0.9, tCov=0.5, fident=0.5),                                    True,  "medium"),
    (_row(qCov=0.9, tCov=0.5, fident=0.0, evalue="1e-7"),                     True,  "medium"),
    (_row(qCov=0.1, tCov=0.1, fident=0.1, evalue="1e-2"),                     True,  "low"),
    # no-structures (prostt5) path
    (_row(qCov=0.9, tCov=0.9, fident=0.0, evalue="1e-3", prostt5_confidence=70),  False, "high"),
    (_row(qCov=0.9, tCov=0.5, fident=0.0, evalue="1e-7", prostt5_confidence=50),  False, "medium"),
    (_row(qCov=0.9, tCov=0.5, fident=0.0, evalue="1e-3", prostt5_confidence=50),  False, "low"),
    (_row(qCov=0.1, tCov=0.1, fident=0.1, evalue="1e-2", prostt5_confidence=10),  False, "low"),
    # boundary cases (mirror the parametrised set above)
    (_row(qCov=0.8,  tCov=0.9,  fident=0.5,  evalue="1e-20"),                          True,  "medium"),
    (_row(qCov=0.81, tCov=0.81, fident=0.31, evalue="1"),                              True,  "high"),
    (_row(qCov=0.81, tCov=0.81, fident=0.3,  evalue="1e-11"),                          True,  "high"),
    (_row(qCov=0.81, tCov=0.81, fident=0.3,  evalue="1e-10"),                          True,  "medium"),
    (_row(qCov=0.9,  tCov=0.9,  fident=0.0,  evalue="1",     prostt5_confidence=61),   False, "high"),
    (_row(qCov=0.9,  tCov=0.9,  fident=0.0,  evalue="1",     prostt5_confidence=60),   False, "low"),
    (_row(qCov=0.9,  tCov=0.5,  fident=0.0,  evalue="1e-7",  prostt5_confidence=45),   False, "medium"),
    (_row(qCov=0.9,  tCov=0.5,  fident=0.0,  evalue="1e-7",  prostt5_confidence=44),   False, "low"),
]


@pytest.mark.parametrize("structures", [True, False])
def test_vectorised_matches_scalar(structures):
    """Apply the polars expression to a stacked DataFrame of all scalar-test
    rows for this ``structures`` value and confirm every row gets the same
    label that the scalar function would assign. Catches any drift between
    the two implementations the moment it appears."""
    rows = [r for (r, s, _exp) in _AGREEMENT_ROWS if s == structures]
    expected = [exp for (_r, s, exp) in _AGREEMENT_ROWS if s == structures]

    pl_df = pl.DataFrame(rows)
    out = pl_df.with_columns(
        assign_annotation_confidence_expr(structures=structures)
        .alias("annotation_confidence")
    )
    actual = out["annotation_confidence"].to_list()
    assert actual == expected, (
        f"vectorised expression diverged from scalar (structures={structures})\n"
        f"  expected: {expected}\n"
        f"  actual:   {actual}"
    )


# ─── zero Foldseek hits through the whole compare pipeline ──────────────────
# Regression for v1.3.0 (phage-annotation.org): a 1.7 kb contig, a 23S rRNA
# fragment, gave no Foldseek hits and phold crashed with
# ``polars.exceptions.NoDataError: empty CSV``. These tests drive everything
# after the Foldseek search (topfunction -> GenBank -> per-CDS TSV -> sub-DB
# TSVs -> function counts) via the --restart path, which reads an existing
# foldseek_results.tsv instead of running Foldseek. No Foldseek binary,
# ProstT5 model or phold database is needed.

TEST_DATA = Path(__file__).parent.parent / "test_data"
# pharokka v1.11.0 (-g phanotate --fast) output for that contig, which is
# tests/test_data/no_foldseek_hits.fasta. Every CDS is No_PHROGs_HMM.
NO_HITS_GBK = TEST_DATA / "no_foldseek_hits_pharokka.gbk"
# 7 of its 17 CDS have a pharokka PHROG, which must survive zero hits.
NC_043029_GBK = TEST_DATA / "NC_043029_pharokka1.4.1.gbk"
PREFIX = "phold"
PROSTT5_CONFIDENCE = 55.5


def _cds_features(gbk: Path) -> dict:
    """{(contig_id, cds_id): CDS feature} for a pharokka GenBank."""
    return {
        (record.id, feature.qualifiers["ID"][0]): feature
        for record in SeqIO.parse(gbk, "genbank")
        for feature in record.features
        if feature.type == "CDS"
    }


def _compare_with_no_foldseek_hits(
    gbk: Path,
    output: Path,
    database: Path,
    structures: bool = False,
    proteins_flag: bool = False,
) -> bool:
    """Run ``subcommand_compare --restart`` on ``gbk`` with an empty
    foldseek_results.tsv, writing the other files --restart reads first."""
    output.mkdir()
    cds = _cds_features(gbk)

    if proteins_flag:
        # same shape as the proteins-compare CLI builds
        gb_dict = {
            "proteins": {
                cds_id: SeqFeature(
                    FeatureLocation(0, len(f.qualifiers["translation"][0])),
                    type="CDS",
                    qualifiers={"ID": cds_id, "translation": f.qualifiers["translation"][0]},
                )
                for (_, cds_id), f in cds.items()
            }
        }
    else:
        gb_dict, _ = get_genbank(str(gbk))

    (output / "foldseek_results.tsv").write_text("")
    if not structures:
        # --restart only rewrites '~PIPE~' in these, so the 3Di content is moot
        aa = "".join(
            f">{cds_id if proteins_flag else f'{contig_id}:{cds_id}'}\n"
            f"{f.qualifiers['translation'][0]}\n"
            for (contig_id, cds_id), f in cds.items()
        )
        (output / f"{PREFIX}_aa.fasta").write_text(aa)
        (output / f"{PREFIX}_3di.fasta").write_text(aa)
        (output / f"{PREFIX}_prostT5_3di_mean_probabilities.csv").write_text(
            "".join(f"{cds_id},{PROSTT5_CONFIDENCE}\n" for _, cds_id in cds)
        )

    return subcommand_compare(
        gb_dict,
        output,
        threads=1,
        evalue=1e-3,
        card_vfdb_evalue=1e-10,
        sensitivity=9.5,
        database=database,
        prefix=PREFIX,
        predictions_dir=output,
        structures=structures,
        structure_dir=output if structures else None,
        logdir=output / "logs",
        filter_structures=False,
        remote_flag=True,
        proteins_flag=proteins_flag,
        fasta_flag=False,
        separate=False,
        max_seqs=1000,
        ultra_sensitive=False,
        extra_foldseek_params=None,
        custom_db=None,
        foldseek_gpu=False,
        restart=True,
    )


def _read_per_cds(output: Path) -> pl.DataFrame:
    # all-Utf8 so empty Foldseek fields come back as null regardless of dtype
    return pl.read_csv(
        output / f"{PREFIX}_per_cds_predictions.tsv",
        separator="\t",
        infer_schema_length=0,
    )


@pytest.mark.parametrize("gbk", [NO_HITS_GBK, NC_043029_GBK], ids=lambda p: p.stem)
def test_compare_no_foldseek_hits_keeps_input_annotations(
    tmp_path, mini_phold_db, loguru_errors, gbk
):
    output = tmp_path / "out"
    assert _compare_with_no_foldseek_hits(gbk, output, mini_phold_db) is True

    # ERROR exits phold, so zero hits must only warn
    assert loguru_errors == []

    cds = _cds_features(gbk)
    per_cds = _read_per_cds(output)
    assert per_cds.height == len(cds)
    for row in per_cds.iter_rows(named=True):
        feature = cds[(row["contig_id"], row["cds_id"])]
        phrog = feature.qualifiers["phrog"][0].replace("No_PHROGs_HMM", "No_PHROG")
        method = "none" if phrog == "No_PHROG" else "pharokka"
        assert row["phrog"] == phrog
        assert row["product"] == feature.qualifiers["product"][0]
        assert row["annotation_method"] == method
        assert row["annotation_confidence"] == method
        assert row["bitscore"] is None
        assert row["tophit_protein"] is None
        assert float(row["prostt5_confidence"]) == PROSTT5_CONFIDENCE

    # GenBank: same CDS with the input products, non-CDS features kept
    in_records = list(SeqIO.parse(gbk, "genbank"))
    out_records = list(SeqIO.parse(output / f"{PREFIX}.gbk", "genbank"))
    assert [r.id for r in out_records] == [r.id for r in in_records]
    for in_rec, out_rec in zip(in_records, out_records):
        assert sorted(f.type for f in out_rec.features) == sorted(f.type for f in in_rec.features)
        out_products = {
            f.qualifiers["ID"][0]: f.qualifiers["product"][0]
            for f in out_rec.features if f.type == "CDS"
        }
        assert out_products == {
            f.qualifiers["ID"][0]: f.qualifiers["product"][0]
            for f in in_rec.features if f.type == "CDS"
        }

    counts = pl.read_csv(output / f"{PREFIX}_all_cds_functions.tsv", separator="\t")
    assert counts.filter(pl.col("Description") == "CDS")["Count"].sum() == len(cds)

    for subdb in ("acr", "vfdb", "card", "netflax", "defensefinder"):
        assert (output / "sub_db_tophits" / f"{subdb}_cds_predictions.tsv").exists()


@pytest.mark.parametrize(
    "structures,proteins_flag",
    [(True, False), (False, True), (True, True)],
    ids=["structures", "proteins", "proteins_structures"],
)
def test_compare_no_foldseek_hits_other_modes(
    tmp_path, mini_phold_db, loguru_errors, structures, proteins_flag
):
    """compare --structures and proteins-compare also finish with zero hits."""
    output = tmp_path / "out"
    assert _compare_with_no_foldseek_hits(
        NO_HITS_GBK, output, mini_phold_db, structures=structures, proteins_flag=proteins_flag
    ) is True
    assert loguru_errors == []

    per_cds = _read_per_cds(output)
    assert sorted(per_cds["cds_id"]) == sorted(cds_id for _, cds_id in _cds_features(NO_HITS_GBK))
    assert per_cds["bitscore"].null_count() == per_cds.height
    assert set(per_cds["product"]) == {"hypothetical protein"}
