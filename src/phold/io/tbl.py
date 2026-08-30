"""
NCBI feature-table (``.tbl``) writer — GenBank submission output.

Mirrors the ``.tbl`` that Pharokka emits (``pharokka.post_processing.create_tbl``)
so a Phold-reannotated genome can be handed to ``table2asn`` for GenBank
submission without a format conversion step (issue #137).

Format (tab-separated, per the NCBI feature-table spec):

    >Feature <contig_id>
    <start>\t<stop>\t<feature_key>
    \t\t\t<qualifier>\t<value>

``start``/``stop`` are 1-based and inclusive, and are written in *transcription*
order: for a minus-strand feature the larger coordinate comes first. That
ordering is how the feature table encodes strand — there is no separate strand
column — so it is the one piece of the format that must not be "tidied".
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import polars as pl
from loguru import logger

from Bio.SeqFeature import SeqFeature

from phold.utils.util import atomic_write_path


# Qualifiers emitted verbatim for a CDS, in the order NCBI's examples use.
# Values come from per_cds_df, which write_genbank() has already normalised to
# 1-based inclusive coordinates (see the issue #75/#77 comment there).
# ``inference`` is not here — it has to be built (see _format_inference).
_CDS_QUALIFIER_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("product", "product"),
    ("function", "function"),
)

# INSDC requires /inference to begin with a category from a controlled
# vocabulary, followed by ":<database>:<accession>". Phold's raw
# annotation_method values ("foldseek", "pharokka", "none") have no category,
# so table2asn rejects every one of them outright with "Qualifier had bad
# value" — one hard error per CDS.
#
# Category choice, measured against table2asn 1.28.1179:
#
#   similar to AA sequence              -> InvalidInferenceValue warning
#   similar to AA sequence:PHROG:1215   -> InvalidInferenceValue warning
#                                          ("unrecognized database")
#   protein motif:PHROG:1215            -> clean
#   alignment:Foldseek:1215             -> clean
#
# The "similar to ... sequence" categories are checked against NCBI's list of
# recognised databases, which PHROG is not on. "protein motif" is not, and is
# also the honest description: a PHROG is a protein group / profile, and Phold
# transfers its annotation from the group its structure matched.
_INFERENCE_CATEGORY = "protein motif:PHROG"

# annotation_method values that mean "no annotation was transferred", for
# which no evidence can honestly be cited.
_NO_INFERENCE_METHODS = frozenset({"none", "", "no_phrog"})


# Phold/Pharokka bookkeeping qualifiers that must never reach a submission
# file: either they are internal identifiers, or table2asn rejects them.
# ``ID``/``locus_tag`` are dropped because NCBI assigns locus tags at
# submission time; ``score``/``source``/``phase`` are pipeline provenance.
_NON_CDS_QUALIFIER_BLOCKLIST = frozenset(
    {
        "ID",
        "locus_tag",
        "score",
        "source",
        "phase",
        "translation",
        "transl_table",
        "phrog",
        "top_hit",
        "isotype",
    }
)

# Pharokka writes the tRNA name to a ``/trna`` qualifier, which is not a valid
# NCBI qualifier — ``product`` is, and it is mandatory for the tRNA feature key.
# ``write_genbank`` renames it in place, but only for records it walks, so do
# the same mapping here rather than depending on that call having happened
# first: a tRNA reaching table2asn without a product is rejected.
_QUALIFIER_ALIASES: Tuple[Tuple[str, str], ...] = (("trna", "product"),)

# NCBI's /anticodon must be a location, e.g. "(pos:115194..115196,aa:Met)".
# Pharokka < v1.8.2 wrote a bare 3-letter code ("CAT"), which table2asn
# rejects with "Qualifier had bad value" — write_genbank already warns about
# that input, so drop the qualifier here rather than emit a known-invalid one.
_ANTICODON_LOCATION = re.compile(r"^\(.*pos:.*\)$", re.IGNORECASE)

# Preferred emission order for the non-CDS qualifiers Pharokka carries. Anything
# not listed keeps dict order after these.
_NON_CDS_QUALIFIER_ORDER: Tuple[str, ...] = (
    "gene",
    "product",
    "inference",
    "ncRNA_class",
    "anticodon",
    "tag_peptide",
    "rpt_family",
    "rpt_type",
    "rpt_unit_range",
    "rpt_unit_seq",
    "db_xref",
    "note",
)


def _format_inference(method: Optional[str], phrog: Optional[str]) -> Optional[str]:
    """Build a valid /inference string, or None to omit the qualifier.

    Cites the PHROG the annotation came from. Without one there is no evidence
    to point at, and a bare category is itself flagged by table2asn, so the
    qualifier is omitted rather than padded with something invented.
    """
    if method is None or str(method).strip().lower() in _NO_INFERENCE_METHODS:
        return None

    if phrog is None:
        return None

    # "No_PHROG" is Phold's sentinel for a CDS with no PHROG assignment.
    text = str(phrog).strip()
    if not text or text.lower() == "no_phrog":
        return None

    return f"{_INFERENCE_CATEGORY}:{text}"


def _is_minus(strand: Optional[Union[int, str]]) -> bool:
    """True for a minus-strand feature.

    Strand reaches this module in two representations: ``per_cds_df`` has
    already been converted to "+"/"-" strings by ``write_genbank``, while the
    non-CDS ``SeqFeature`` objects still carry BioPython's ±1 ints. Accept
    both — reading one as the other silently mis-orients every reverse-strand
    feature in the file.
    """
    if strand is None:
        return False
    if isinstance(strand, str):
        return strand.strip() in {"-", "-1"}
    try:
        return int(strand) == -1
    except (TypeError, ValueError):
        return False


def _orient(
    start: int, end: int, strand: Optional[Union[int, str]]
) -> Tuple[int, int]:
    """Return (start, stop) in transcription order.

    A minus-strand feature is written high-coordinate first; that swap *is*
    the strand annotation in a feature table.
    """
    if _is_minus(strand):
        return end, start
    return start, end


def _split_partial(
    partial: Optional[str], strand: Optional[Union[int, str]]
) -> Tuple[bool, bool]:
    """Map a Prodigal ``partial`` flag onto (5'-incomplete, 3'-incomplete).

    ``partial`` is two digits in *genomic* orientation — left edge then right
    edge — regardless of strand (verified against pyrodigal's own GFF writer:
    reverse-complementing a sequence turns a left-edge ``10`` gene into a
    right-edge ``01`` gene). So the 5' end is the left digit on the plus
    strand and the right digit on the minus strand.
    """
    if partial is None:
        return False, False

    text = str(partial).strip()
    if len(text) != 2 or set(text) - {"0", "1"}:
        return False, False

    left, right = text[0] == "1", text[1] == "1"
    if _is_minus(strand):
        return right, left
    return left, right


def _cds_coordinates(
    five_end: int,
    three_end: int,
    strand: Optional[Union[int, str]],
    partial: Optional[str],
    contig_length: Optional[int],
) -> Tuple[str, str, Optional[int]]:
    """Render a CDS's coordinate pair, with partial markers and codon_start.

    ``five_end``/``three_end`` are taken straight from ``per_cds_df``, which is
    **already in transcription order**: ``write_genbank`` assigns
    ``start = location.end`` for a minus-strand CDS, so start > end there. They
    must not be re-oriented — doing so writes every reverse-strand CDS
    backwards. (The non-CDS path is different: those come from BioPython
    locations, which really are genomic min/max, hence ``_orient``.)

    NCBI marks an incomplete 5' end with ``<`` on the first coordinate and an
    incomplete 3' end with ``>`` on the second. When the 5' end is incomplete
    but does not sit flush against the contig edge, the feature is extended to
    the edge and ``codon_start`` carries the reading-frame offset.

    Returns ``(start_text, stop_text, codon_start)``.
    """
    five_partial, three_partial = _split_partial(partial, strand)
    start_text, stop_text = str(five_end), str(three_end)
    codon_start: Optional[int] = None

    if three_partial:
        # NCBI expects a 3'-partial feature to run to the sequence edge. Gene
        # callers stop at the last whole codon, leaving 1-2 trailing bases, and
        # table2asn then warns "3' partial is not at end of sequence". Extend
        # to the edge, mirroring the 5' handling below.
        edge = 1 if _is_minus(strand) else contig_length
        if edge is not None:
            three_end = edge
        stop_text = f">{three_end}"

    if five_partial:
        start_text = f"<{five_end}"

        # Distance between the 5' end and the contig edge it ran off. Anything
        # non-zero means the annotation must be pushed out to the edge, with
        # the leftover bases declared as the frame offset. The minus-strand 5'
        # end is the high coordinate, so it runs off the far end of the contig.
        if _is_minus(strand):
            if contig_length is not None and contig_length - five_end > 0:
                codon_start = contig_length - five_end + 1
                start_text = f"<{contig_length}"
        elif five_end > 1:
            codon_start = five_end
            start_text = "<1"

        # codon_start is a frame offset, so only 1/2/3 are meaningful. Warn
        # rather than abort: this runs at the very end of a long pipeline, and
        # a warning the user can act on beats losing the whole run's output.
        if codon_start is not None and codon_start > 3:
            logger.warning(
                f"codon_start of {codon_start} for the partial CDS at "
                f"{five_end}..{three_end} is greater than 3, which NCBI will "
                "reject. Please raise an issue on GitHub with your genome."
            )

    return start_text, stop_text, codon_start


def _write_cds_features(
    handle, contig_df: pl.DataFrame, contig_length: Optional[int]
) -> int:
    """Write every CDS row of one contig. Returns the number written."""
    written = 0
    for row in contig_df.iter_rows(named=True):
        start_text, stop_text, codon_start = _cds_coordinates(
            int(row["start"]),
            int(row["end"]),
            row.get("strand"),
            row.get("partial"),
            contig_length,
        )
        handle.write(f"{start_text}\t{stop_text}\tCDS\n")

        for column, qualifier in _CDS_QUALIFIER_COLUMNS:
            value = row.get(column)
            # A null product/function would emit a bare qualifier line that
            # table2asn treats as a parse error, so skip rather than write "".
            if value is None or str(value).strip() == "":
                continue
            handle.write(f"\t\t\t{qualifier}\t{value}\n")

        inference = _format_inference(row.get("annotation_method"), row.get("phrog"))
        if inference is not None:
            handle.write(f"\t\t\tinference\t{inference}\n")

        transl_table = row.get("transl_table")
        if transl_table is not None and str(transl_table).strip() != "":
            handle.write(f"\t\t\ttransl_table\t{transl_table}\n")

        # After transl_table, matching Pharokka's field order.
        if codon_start is not None:
            handle.write(f"\t\t\tcodon_start\t{codon_start}\n")
        written += 1
    return written


def _non_cds_sort_key(item: Tuple[str, SeqFeature]) -> int:
    """Sort non-CDS features by start coordinate, tolerating odd locations."""
    _, feature = item
    try:
        return int(feature.location.start)
    except (AttributeError, TypeError):
        return 0


def _write_non_cds_features(
    handle, features: Dict[str, SeqFeature]
) -> int:
    """Write tRNA / tmRNA / CRISPR / ncRNA features carried over from the input.

    Phold does not predict these — they are passed through from the Pharokka
    (or NCBI/Bakta) input GenBank — but a submission-ready feature table has to
    carry them, otherwise the .tbl silently describes fewer features than the
    .gbk Phold wrote alongside it.
    """
    written = 0
    for _, feature in sorted(features.items(), key=_non_cds_sort_key):
        location = getattr(feature, "location", None)
        if location is None:
            continue

        # BioPython locations are 0-based half-open; the feature table is
        # 1-based inclusive, so the start needs +1 and the end is already
        # the inclusive last base.
        try:
            start = int(location.start) + 1
            end = int(location.end)
        except (TypeError, ValueError):
            logger.warning(
                f"Skipping feature with an unparseable location in the .tbl: {feature.type}"
            )
            continue

        start, stop = _orient(start, end, location.strand)
        handle.write(f"{start}\t{stop}\t{feature.type}\n")

        # Copy: the same SeqFeature objects are shared with the GenBank writer,
        # so the alias rewrite below must not mutate the caller's features.
        qualifiers = dict(feature.qualifiers or {})
        for source_key, target_key in _QUALIFIER_ALIASES:
            aliased = qualifiers.pop(source_key, None)
            if aliased is not None and target_key not in qualifiers:
                qualifiers[target_key] = aliased

        ordered = [k for k in _NON_CDS_QUALIFIER_ORDER if k in qualifiers]
        ordered += [k for k in qualifiers if k not in _NON_CDS_QUALIFIER_ORDER]

        for key in ordered:
            if key in _NON_CDS_QUALIFIER_BLOCKLIST:
                continue
            value = qualifiers[key]
            # BioPython stores qualifiers as lists when parsed from a GenBank
            # file but as bare strings when Phold built the feature itself.
            if isinstance(value, (list, tuple)):
                value = value[0] if value else None
            if value is None or str(value).strip() == "":
                continue
            if key == "anticodon" and not _ANTICODON_LOCATION.match(str(value).strip()):
                logger.warning(
                    f"Dropping non-location anticodon '{value}' from the .tbl — "
                    "re-run Pharokka >= v1.8.2 to submit this tRNA to GenBank."
                )
                continue
            handle.write(f"\t\t\t{key}\t{value}\n")
        written += 1
    return written


def write_tbl(
    per_cds_df: pl.DataFrame,
    non_cds_dict: Dict[str, Dict[str, SeqFeature]],
    contig_ids: List[str],
    prefix: str,
    output: Path,
    contig_lengths: Optional[Dict[str, int]] = None,
) -> Path:
    """Write the NCBI feature table for every contig.

    Args:
        per_cds_df: The per-CDS dataframe returned by ``write_genbank`` —
            already 1-based inclusive, with contig_id / start / end / strand /
            product / function / annotation_method / transl_table, and
            ``partial`` when the input carried Prodigal partial flags.
        non_cds_dict: ``{contig_id: {feature_id: SeqFeature}}`` for the
            tRNA / tmRNA / CRISPR features passed through from the input.
        contig_ids: Contig order to emit, so the .tbl matches the .gbk.
        prefix: Output filename prefix.
        output: Output directory.
        contig_lengths: ``{contig_id: length}``, needed to place a minus-strand
            5'-partial CDS against the contig edge. Without it those CDS still
            get their ``<`` marker, just no codon_start.

    Returns:
        Path to the written .tbl.

    Note:
        Every contig gets a ``>Feature`` header even when it has no features,
        because table2asn matches records by that header — omitting an empty
        contig makes it look absent from the submission rather than empty.
    """
    out_path: Path = Path(output) / f"{prefix}.tbl"

    has_cds = per_cds_df.height > 0 and "contig_id" in per_cds_df.columns
    total = 0

    with atomic_write_path(out_path) as tmp, open(tmp, "w") as f:
        for contig_id in contig_ids:
            f.write(f">Feature {contig_id}\n")

            if has_cds:
                contig_df = per_cds_df.filter(pl.col("contig_id") == contig_id)
                total += _write_cds_features(
                    f, contig_df, (contig_lengths or {}).get(contig_id)
                )

            total += _write_non_cds_features(f, non_cds_dict.get(contig_id, {}))

    logger.info(f"Wrote {total} features across {len(contig_ids)} contig(s) to {out_path}")
    return out_path
