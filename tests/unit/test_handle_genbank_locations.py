"""
Function-level tests for `write_genbank`'s handling of unparseable feature locations.

BioPython sets `SeqFeature.location` to None when it cannot parse a location, and warns
instead of raising:

    BiopythonParserWarning: negative starting position in feature location '-63..456';
    setting feature location to None

Aragorn emits such coordinates for a tmRNA that runs off the start of a contig, and Pharokka
writes them into the Genbank unchanged, so phold receives them through ordinary input rather
than a corrupt file. `write_genbank` sorted features with `key=lambda x: x.location.start`,
which raised `AttributeError: 'NoneType' object has no attribute 'start'` at the final write
step -- discarding an entire run's ProstT5 and Foldseek predictions for one bad tRNA feature.

Found on 6 of ~540 PhageScope chunks (~1.1%), each time a tmRNA.
"""

from __future__ import annotations

import warnings

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from phold.io.handle_genbank import write_genbank


def _cds(start: int, end: int, name: str) -> SeqFeature:
    feature = SeqFeature(FeatureLocation(start, end, strand=1), type="CDS")
    feature.qualifiers = {
        "ID": [name],
        "phrog": ["1"],
        "function": ["tail"],
        "product": ["tail protein"],
        "translation": ["MA"],
        "transl_table": ["11"],
    }
    return feature


def _unparseable(name: str) -> SeqFeature:
    """A feature as BioPython hands it back when the location could not be parsed."""
    feature = SeqFeature(None, type="tmRNA")
    feature.qualifiers = {"locus_tag": [name]}
    return feature


@pytest.fixture
def record() -> dict[str, SeqRecord]:
    seq = SeqRecord(Seq("ATG" * 200), id="contig_1", name="contig_1", description="")
    seq.annotations = {"molecule_type": "DNA"}
    return {"contig_1": seq}


def test_unparseable_location_does_not_abort_the_write(tmp_path, record):
    """One None location must not cost the whole contig."""
    updated = {"contig_1": {"cds_1": _cds(0, 60, "contig_1_CDS_0001")}}
    non_cds = {"contig_1": {"tmrna_1": _unparseable("contig_1_tmRNA_0001")}}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        per_cds = write_genbank(
            updated_cds_dict=updated,
            non_cds_dict=non_cds,
            source_dict={"contig_1": {"contig_1_CDS_0001": "foldseek"}},
            prefix="test",
            gb_dict=record,
            output=tmp_path,
            proteins_flag=False,
            separate=False,
            fasta_flag=False,
        )

    # The CDS survives; only the feature BioPython could not place is dropped.
    assert len(per_cds) == 1
    assert (tmp_path / "test.gbk").exists()


def test_features_with_locations_are_still_sorted(tmp_path, record):
    """Dropping the unplaceable feature must not disturb ordering of the rest."""
    updated = {
        "contig_1": {
            "cds_2": _cds(120, 180, "contig_1_CDS_0002"),
            "cds_1": _cds(0, 60, "contig_1_CDS_0001"),
        }
    }
    non_cds = {"contig_1": {"tmrna_1": _unparseable("contig_1_tmRNA_0001")}}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        per_cds = write_genbank(
            updated_cds_dict=updated,
            non_cds_dict=non_cds,
            source_dict={
                "contig_1": {
                    "contig_1_CDS_0001": "foldseek",
                    "contig_1_CDS_0002": "foldseek",
                }
            },
            prefix="test",
            gb_dict=record,
            output=tmp_path,
            proteins_flag=False,
            separate=False,
            fasta_flag=False,
        )

    assert per_cds["start"].to_list() == sorted(per_cds["start"].to_list())
