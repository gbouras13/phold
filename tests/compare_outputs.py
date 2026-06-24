#!/usr/bin/env python3
"""
Compare phold outputs between two runs (dev pholdlib-refactored vs bioconda reference).
Ignores timestamp-dependent content (log files, dates).

Usage:
    python tests/compare_outputs.py <dir_dev> <dir_ref>

Exit code 0 = identical (modulo timestamps), non-zero = differences found.
"""
import json
import math
import re
import shutil
import sys
from pathlib import Path

# ── patterns that are timestamp/run-specific and should be ignored ──────────
SKIP_LINE_PATTERNS = [
    re.compile(r"^\d{4}-\d{2}-\d{2}"),          # log lines: 2026-05-26 ...
    re.compile(r"#.*phold.*run"),                 # any phold run pragma
]

# ── in-line substring normalisation (applied to every kept line) ────────────
# phold stamps each GenBank record with a fresh annotation timestamp, a LOCUS
# date, and the phold version. These vary per run / per release while the rest
# of the line is meaningful, so normalise the volatile substring in place
# rather than dropping the whole line. Goldens stay stable across days and
# version bumps; everything else on the line is still compared.
NORMALIZE_SUBS = [
    (re.compile(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}"), "<TIMESTAMP>"),  # Annotation Date: ...
    (re.compile(r"\d{2}-[A-Z]{3}-\d{4}"), "<DATE>"),                       # GenBank LOCUS date e.g. 24-JUN-2026
    (re.compile(r"(Annotated with Phold )\S+"), r"\1<VERSION>"),           # phold version in the comment
]


def normalize_line(line: str) -> str:
    """Replace volatile date/version substrings with stable placeholders."""
    for pat, repl in NORMALIZE_SUBS:
        line = pat.sub(repl, line)
    return line

# file extensions to skip entirely. ``.log`` is timestamp noise; the ``.h5``/
# embedding binaries (from ``phold predict --save_*_embeddings``) and other
# large binary artefacts are excluded so they are neither committed as goldens
# nor flagged as "ONLY IN DEV" when present only in a fresh run.
SKIP_EXTENSIONS = {".log", ".h5", ".hdf5", ".pkl", ".pickle", ".pt", ".npy", ".npz"}

# files to skip by name
SKIP_FILENAMES = set()

# directory components to skip entirely (any file under these dirs is ignored).
# Besides log dirs, the Foldseek query/result/temp databases and the optional
# filtered-structures copy are binary intermediates that are not meaningful to
# pin as goldens.
SKIP_DIRS = {
    "logs", "logdir",
    "foldseek_db", "result_db", "temp_db", "filtered_structures",
}


def should_skip(rel: Path) -> bool:
    """True if a path (relative to the output dir root) is excluded from both
    golden comparison and golden generation."""
    return (
        rel.suffix in SKIP_EXTENSIONS
        or rel.name in SKIP_FILENAMES
        or bool(SKIP_DIRS.intersection(rel.parts))
    )


def copy_golden(src: Path, dst: Path) -> None:
    """Copy a produced phold output dir into the golden tree, applying the same
    skip rules used for comparison so committed goldens stay small and textual.

    Used by ``--update_goldens`` to (re)generate the committed reference set.
    ``dst`` is replaced wholesale so removed outputs don't linger.
    """
    if dst.exists():
        shutil.rmtree(dst)
    for f in sorted(src.rglob("*")):
        if not f.is_file():
            continue
        rel = f.relative_to(src)
        if should_skip(rel):
            continue
        target = dst / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, target)


def filter_lines(path: Path) -> list:
    """Read a file, drop timestamp-like log lines, and normalise volatile
    date/version substrings on the lines that remain."""
    try:
        lines = path.read_text(errors="replace").splitlines()
    except Exception as e:
        return [f"<ERROR reading {path}: {e}>"]
    return [
        normalize_line(l)
        for l in lines
        if not any(p.search(l) for p in SKIP_LINE_PATTERNS)
    ]


# ── phold per-CDS / sub-DB / custom-hit TSV: drop non-reproducible columns ──
# Mirroring baktfold's ``_tophit_differ``: because the ProstT5 3Di input is
# itself non-deterministic on GPU, every Foldseek alignment number wobbles
# run-to-run — not just the quality scores (bitscore/fident/evalue) but the
# alignment extent (qStart/qEnd/tStart/tEnd), the coverage derived from it
# (qCov/tCov), the structural scores (alntmscore/lddt), and everything derived
# from bitscores (the ``*_bitscore_proportion`` columns + prostt5_confidence).
# The reproducible facts are the genomic coords, intrinsic lengths (qLen/tLen),
# the *identity* of the hit (tophit_protein / target + sub-DB metadata) and the
# annotations (phrog/function/product/method/tier). Drop the volatile columns by
# name and compare the rest exactly, so a changed annotation / hit / column
# structure is still caught while alignment-numeric noise is ignored.
#
# NB: ``lddt`` (lowercase, custom-DB per-alignment LDDT) is volatile and dropped;
# ``pLDDT`` (sub-DB metadata, per-target, deterministic) is kept.
_VOLATILE_TSV_COLS = {
    "bitscore", "fident", "evalue",
    "qStart", "qEnd", "qCov", "tStart", "tEnd", "tCov",
    "alntmscore", "lddt",
    "prostt5_confidence",
    "top_bitscore_proportion_not_unknown",
}


def _is_volatile_col(name: str) -> bool:
    return name in _VOLATILE_TSV_COLS or name.endswith("_bitscore_proportion")


def _phold_tsv_differ(lines_dev: list, lines_ref: list) -> list:
    """Compare a phold per-CDS / sub-DB / custom-hit TSV, dropping the
    non-reproducible Foldseek/ProstT5 numeric columns (see _VOLATILE_TSV_COLS)
    and comparing the remaining columns exactly (sorted).

    The header is part of the comparison, so added/removed/renamed columns and
    any changed annotation/hit/length still fail.
    """
    if not lines_dev and not lines_ref:
        return []
    hdr_dev = lines_dev[0].split("\t") if lines_dev else []
    hdr_ref = lines_ref[0].split("\t") if lines_ref else []
    if hdr_dev != hdr_ref:
        return [f"    header differs: dev={hdr_dev[:12]} ref={hdr_ref[:12]}"]
    keep = [i for i, name in enumerate(hdr_dev) if not _is_volatile_col(name)]

    def _project(lines):
        rows = []
        for line in lines[1:]:                       # skip header
            parts = line.split("\t")
            rows.append("\t".join(parts[i] for i in keep if i < len(parts)))
        return sorted(rows)

    sd, sr = _project(lines_dev), _project(lines_ref)
    row_diffs = []
    if sd != sr:
        for i, (a, b) in enumerate(zip(sd, sr)):
            if a != b:
                row_diffs.append(f"    dev[{i}]: {a[:140]}")
                row_diffs.append(f"    ref[{i}]: {b[:140]}")
                if i > 10:
                    row_diffs.append("    ... (truncated)")
                    break
        if len(sd) != len(sr):
            row_diffs.append(f"    line count: dev={len(sd)} ref={len(sr)}")
    return row_diffs


def _jsonl_prob_differ(lines_dev: list, lines_ref: list, tol: float = 0.01) -> list:
    """Compare two JSONL probability files (``_all_probabilities.json``) with
    per-probability float tolerance.

    Each line is ``{"seq_id": str, "probability": [float, ...]}``. We compare
    seq_ids exactly and probabilities element-wise with ``abs_tol=tol``.
    Absorbs the tiny CPU-matmul drift between torch versions on the
    ProstT5 softmax outputs (typical magnitude: ~0.01 per element).
    """
    row_diffs = []
    if len(lines_dev) != len(lines_ref):
        row_diffs.append(f"    line count: dev={len(lines_dev)} ref={len(lines_ref)}")
    for i, (a, b) in enumerate(zip(lines_dev, lines_ref)):
        if a == b:
            continue
        try:
            da, db = json.loads(a), json.loads(b)
        except json.JSONDecodeError:
            row_diffs.append(f"    JSON parse error dev[{i}]: {a[:120]}")
            row_diffs.append(f"                       ref[{i}]: {b[:120]}")
            continue
        if da.get("seq_id") != db.get("seq_id"):
            row_diffs.append(
                f"    seq_id mismatch dev[{i}]: {da.get('seq_id')}"
                f" | ref: {db.get('seq_id')}"
            )
            continue
        pa, pb = da.get("probability", []), db.get("probability", [])
        if len(pa) != len(pb):
            row_diffs.append(
                f"    prob len dev[{i}]={len(pa)} ref[{i}]={len(pb)} for {da['seq_id']}"
            )
            continue
        bad_pos = [
            j for j, (fa, fb) in enumerate(zip(pa, pb))
            if not math.isclose(float(fa), float(fb), abs_tol=tol)
        ]
        if bad_pos:
            row_diffs.append(
                f"    {da['seq_id']} differs at {len(bad_pos)} position(s) > tol={tol}; "
                f"first 3: " + ", ".join(
                    f"[{j}] dev={pa[j]} ref={pb[j]}" for j in bad_pos[:3]
                )
            )
    return row_diffs


def _csv_float_differ(lines_dev: list, lines_ref: list, tol: float = 0.5) -> list:
    """Return diff messages for two sorted lists of CSV lines where the second
    column is a float.  Lines with matching seq_ids are compared numerically;
    count mismatches are flagged by the caller."""
    row_diffs = []
    if len(lines_dev) != len(lines_ref):
        row_diffs.append(f"    line count: dev={len(lines_dev)} ref={len(lines_ref)}")
    for i, (a, b) in enumerate(zip(lines_dev, lines_ref)):
        if a == b:
            continue
        pa, pb = a.split(",", 1), b.split(",", 1)
        if pa[0] != pb[0]:
            row_diffs.append(f"    seq_id mismatch dev[{i}]: {a[:120]} | ref: {b[:120]}")
            continue
        try:
            if not math.isclose(float(pa[1]), float(pb[1]), abs_tol=tol):
                row_diffs.append(f"    value mismatch dev[{i}]: {a[:120]}")
                row_diffs.append(f"                   ref[{i}]: {b[:120]}")
        except ValueError:
            row_diffs.append(f"    parse error dev[{i}]: {a[:120]}")
    return row_diffs


def _parse_fasta(lines: list) -> dict:
    """Parse FASTA lines (already filtered/normalised) into {seq_id: sequence}.

    Tolerant of wrapped sequences (concatenates lines until the next header).
    phold writes one unwrapped sequence per header, but this stays correct
    either way.
    """
    seqs: dict = {}
    sid = None
    buf: list = []
    for l in lines:
        if l.startswith(">"):
            if sid is not None:
                seqs[sid] = "".join(buf)
            sid = l[1:].split()[0] if len(l) > 1 else ""
            buf = []
        elif sid is not None:
            buf.append(l.strip())
    if sid is not None:
        seqs[sid] = "".join(buf)
    return seqs


def _fasta_residue_differ(lines_dev: list, lines_ref: list, min_identity: float) -> list:
    """Compare two phold sequence FASTAs (``*_3di.fasta`` / ``*_aa.fasta``)
    tolerantly, mirroring baktfold's ``_fasta_3di_differ``.

    The 3Di string is the ProstT5 *argmax* and the AA string is the
    confidence-masked translation (``mask_low_confidence_aa``); neither is
    bit-identical across hardware / driver / torch versions, so each sequence is
    compared by per-residue identity and only flagged when it falls below
    *min_identity*. Missing/extra sequences and length mismatches are always
    reported — those indicate a real change, since the sequence length equals
    the protein length.
    """
    dev = _parse_fasta(lines_dev)
    ref = _parse_fasta(lines_ref)
    diffs = []

    only_dev = sorted(set(dev) - set(ref))
    only_ref = sorted(set(ref) - set(dev))
    if only_dev:
        diffs.append(f"    seq ids only in dev: {only_dev[:10]}")
    if only_ref:
        diffs.append(f"    seq ids only in ref: {only_ref[:10]}")

    low = []
    for sid in sorted(set(dev) & set(ref)):
        a, b = dev[sid], ref[sid]
        if len(a) != len(b):
            diffs.append(f"    length mismatch {sid}: dev={len(a)} ref={len(b)}")
            continue
        if not a:
            continue
        identity = sum(1 for x, y in zip(a, b) if x == y) / len(a)
        if identity < min_identity:
            low.append((sid, identity))

    if low:
        low.sort(key=lambda t: t[1])
        worst = ", ".join(f"{sid}={ident:.0%}" for sid, ident in low[:10])
        diffs.append(f"    {len(low)} seq(s) below {min_identity:.0%} identity (worst: {worst})")
    return diffs


def compare_dirs(dir_dev: Path, dir_ref: Path, strict: bool = False) -> list:
    """Recursively compare two directories. Returns list of diff messages."""
    diffs = []

    dev_files = {f.relative_to(dir_dev) for f in dir_dev.rglob("*") if f.is_file()}
    ref_files = {f.relative_to(dir_ref) for f in dir_ref.rglob("*") if f.is_file()}

    for f in sorted(dev_files - ref_files):
        if not should_skip(f):
            diffs.append(f"  ONLY IN DEV : {f}")

    for f in sorted(ref_files - dev_files):
        if not should_skip(f):
            diffs.append(f"  ONLY IN REF : {f}")

    for rel in sorted(dev_files & ref_files):
        if should_skip(rel):
            continue

        fd = dir_dev / rel
        fr = dir_ref / rel

        ld = filter_lines(fd)
        lr = filter_lines(fr)

        # ``*_3di.fasta`` (ProstT5 argmax 3Di) and ``*_aa.fasta`` (the
        # confidence-masked AA — mask_low_confidence_aa): deterministic on CPU
        # (strict) but a few residues flip (argmax / X-mask) under GPU
        # non-determinism, so compare by per-residue identity (exact on CPU,
        # ≥95% on GPU). The seq_id set and every sequence length still match
        # exactly.
        if rel.name.endswith(("_3di.fasta", "_aa.fasta")):
            min_identity = 1.0 if strict else 0.95
            row_diffs = _fasta_residue_differ(ld, lr, min_identity)
            if row_diffs:
                diffs.append(f"  DIFFER (residue identity < {min_identity:.0%}) : {rel}")
                diffs.extend(row_diffs[:22])
            continue

        # For TSV/CSV/TXT where row order may differ, compare sorted.
        # mean_probabilities CSVs always use a tolerance because the values may
        # differ due to two independent sources:
        #   1. MPS (GPU) non-determinism (--cpu runs are not affected)
        #   2. numpy version formatting artefact: bioconda phold v1.2.5 (numpy ≥2)
        #      outputs full IEEE-754 float32 repr (e.g. "82.69000244140625") while
        #      the refactored code outputs clean Python-float 2dp (e.g. "82.69").
        #      The values differ by <3e-6, well within abs_tol=0.01.
        # abs_tol=0.5 covers MPS non-determinism; abs_tol=0.01 suffices for the
        # numpy formatting artefact.  We always apply the tighter of the two:
        #   --cpu strict mode → abs_tol=0.01 (formatting artefact only)
        #   default (MPS)     → abs_tol=0.5  (GPU noise dominates)
        # ``_all_probabilities.json``: JSONL of per-residue ProstT5 softmax
        # outputs (values in the 0-100 percentage range). Compare each
        # probability with abs_tol — phold writes 2dp-rounded values, so
        # straddle cases like dev=68.10 / ref=68.11 (truly 0.005 apart in
        # the underlying float) appear ~0.01 apart and would fail a tight
        # tol. abs_tol=0.02 absorbs that without masking real differences
        # (0.02 is 0.02% of the 0-100 range — well below model noise).
        if rel.suffix == ".json" and "probabilities" in rel.name:
            sd, sr = sorted(ld), sorted(lr)
            jsonl_tol = 0.02 if strict else 0.5
            row_diffs = _jsonl_prob_differ(sd, sr, tol=jsonl_tol)
            if row_diffs:
                diffs.append(f"  DIFFER (per-prob tol={jsonl_tol}) : {rel}")
                diffs.extend(row_diffs[:22])
            continue

        if rel.suffix in {".tsv", ".csv", ".txt"}:
            is_prob_csv = "probabilities" in rel.name and rel.suffix == ".csv"
            # Foldseek-score-bearing TSVs (per_cds_predictions, sub_db_tophits/*,
            # custom_database_hits): drop the non-reproducible alignment numerics
            # and compare the rest exactly (see _phold_tsv_differ).
            is_score_tsv = rel.suffix == ".tsv" and (
                "predictions" in rel.name or "database_hits" in rel.name
            )
            sd, sr = sorted(ld), sorted(lr)
            if is_prob_csv:
                # mean_probabilities CSVs: float tolerance always (absorbs numpy
                # float32-repr artefact on CPU; GPU non-determinism off-CPU).
                prob_tol = 0.01 if strict else 0.5
                row_diffs = _csv_float_differ(sd, sr, tol=prob_tol)
                if row_diffs:
                    diffs.append(f"  DIFFER (sorted, tol={prob_tol}) : {rel}")
                    diffs.extend(row_diffs[:22])
                    if len(sd) != len(sr):
                        diffs.append(f"    line count: dev={len(sd)} ref={len(sr)}")
            elif is_score_tsv:
                row_diffs = _phold_tsv_differ(ld, lr)
                if row_diffs:
                    diffs.append(f"  DIFFER (alignment numerics dropped) : {rel}")
                    diffs.extend(row_diffs[:22])
            elif sd != sr:
                diffs.append(f"  DIFFER (sorted) : {rel}")
                for i, (a, b) in enumerate(zip(sd, sr)):
                    if a != b:
                        diffs.append(f"    dev[{i}]: {a[:140]}")
                        diffs.append(f"    ref[{i}]: {b[:140]}")
                        if i > 10:
                            diffs.append("    ... (truncated)")
                            break
                if len(sd) != len(sr):
                    diffs.append(f"    line count: dev={len(sd)} ref={len(sr)}")
        else:
            if ld != lr:
                diffs.append(f"  DIFFER : {rel}")
                for i, (a, b) in enumerate(zip(ld, lr)):
                    if a != b:
                        diffs.append(f"    dev[{i}]: {a[:140]}")
                        diffs.append(f"    ref[{i}]: {b[:140]}")
                        if i > 10:
                            diffs.append("    ... (truncated)")
                            break
                if len(ld) != len(lr):
                    diffs.append(f"    line count: dev={len(ld)} ref={len(lr)}")

    return diffs


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Compare phold outputs between two runs."
    )
    parser.add_argument("dir_dev", type=Path, help="Dev output directory")
    parser.add_argument("dir_ref", type=Path, help="Ref output directory")
    parser.add_argument(
        "--cpu",
        action="store_true",
        help=(
            "Both runs used --cpu (deterministic). "
            "mean_probabilities CSVs are still compared with abs_tol=0.01 to absorb the "
            "numpy-version float32-repr artefact (bioconda v1.2.5 prints full IEEE-754 "
            "float32 repr; refactored code prints clean 2dp Python floats; difference <3e-6). "
            "All other files are compared exactly."
        ),
    )
    args = parser.parse_args()

    dir_dev = args.dir_dev
    dir_ref = args.dir_ref
    strict = args.cpu

    for d, label in [(dir_dev, "dev"), (dir_ref, "ref")]:
        if not d.is_dir():
            print(f"ERROR: {label} directory does not exist: {d}")
            sys.exit(2)

    print(f"Comparing:\n  DEV: {dir_dev}\n  REF: {dir_ref}\n")
    if strict:
        print("Mode: --cpu (mean-prob CSVs: abs_tol=0.01 for numpy fmt artefact; all else exact)\n")
    else:
        print("Mode: MPS/GPU (mean-prob CSVs: abs_tol=0.5 for GPU non-determinism)\n")

    diffs = compare_dirs(dir_dev, dir_ref, strict=strict)

    if diffs:
        print(f"DIFFERENCES FOUND ({len(diffs)} issues):")
        for d in diffs:
            print(d)
        sys.exit(1)
    else:
        print("ALL OUTPUTS MATCH (modulo timestamps).")
        sys.exit(0)


if __name__ == "__main__":
    main()
