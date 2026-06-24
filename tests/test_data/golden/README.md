# Golden output regression references

This directory holds **committed golden reference outputs** for phold's main
subcommands. Each `tests/test_integration.py` test that produces output also
diffs its result against `golden/<case>/` (via `tests/compare_outputs.py`), so a
future change that silently alters a GenBank annotation, a TSV column, a
confidence tier, or a ProstT5 probability is caught as a test failure.

The golden check is embedded in the integration test itself and runs after
`test_install` has provisioned `tests/test_data/phold_db`, so the database
already exists by the time any golden comparison runs.

## What is pinned

Only the small, meaningful **text** outputs per case — e.g.:

- `<prefix>.gbk` (GenBank)
- `<prefix>_per_cds_predictions.tsv`
- `<prefix>_all_cds_functions.tsv`
- `sub_db_tophits/*.tsv`
- `<prefix>_3di.fasta`, `<prefix>_aa.fasta`
- `<prefix>_prostT5_3di_mean_probabilities.csv`, `*_all_probabilities.json`
- custom-database hit TSVs

Binary / intermediate artefacts are **excluded** by `compare_outputs.py` skip
rules: `logs/`, the Foldseek `foldseek_db`/`result_db`/`temp_db` databases,
`filtered_structures/`, and embedding binaries (`.h5`, …).

## How comparison tolerates GPU non-determinism

ProstT5 is non-deterministic on GPU. That drift cascades into every downstream
output: the argmax 3Di flips a few residues, the confidence-masked AA flips a
few `X`s, and — because the 3Di query changes — the Foldseek alignment numbers
(bitscore/fident/evalue/coords/coverage/structural scores) all wobble too. The
*annotations* (phrog/function/product/method/tier, hit identity, lengths,
counts) stay stable. So, mirroring baktfold's `compare_outputs.py`, the engine
**drops the non-reproducible values and compares the reproducible facts**
(`compare_outputs.compare_dirs(..., strict = not gpu_available)`):

- Per-run timestamps (GenBank `Annotation Date`, `LOCUS` date) and the phold
  version string are normalised in place before comparison.
- `*_3di.fasta` / `*_aa.fasta` are compared by **per-residue identity** — exact
  on CPU (`strict`), ≥95 % on GPU. The seq-id set and every sequence length
  must still match exactly.
- `*_per_cds_predictions.tsv`, `sub_db_tophits/*.tsv` and
  `*_custom_database_hits.tsv`: the volatile Foldseek/ProstT5 numeric columns
  are **dropped by name** (`bitscore`, `fident`, `evalue`, `qStart/qEnd/qCov`,
  `tStart/tEnd/tCov`, `alntmscore`, `lddt`, `prostt5_confidence`, and every
  `*_bitscore_proportion`); the remaining columns — genomic coords, `qLen/tLen`,
  the hit identity (`tophit_protein`/`target` + sub-DB metadata incl. per-target
  `pLDDT`), and the annotations — are compared exactly. A changed annotation,
  hit, length, row count, or column structure still fails.
- ProstT5 probability files (`*_mean_probabilities.csv`, `*_all_probabilities.json`)
  are compared with a float tolerance (tight on CPU, `abs_tol=0.5` on GPU).
- Everything else (`.gbk`, `_all_cds_functions.tsv`, …) is compared exactly,
  modulo the date/version normalisation above.

> **Note:** `annotation_confidence` (high/medium/low) is kept and compared
> exactly. It derives from `prostt5_confidence` thresholds, so a borderline
> value could in principle flip tiers under GPU noise. It was stable across the
> bootstrap run; if it ever flakes, add `annotation_confidence` to
> `_VOLATILE_TSV_COLS` in `tests/compare_outputs.py`.

## Regenerating (after an *intended* output change)

```bash
# after test_install has populated tests/test_data/phold_db
pytest tests/test_integration.py --update_goldens          # all cases
pytest tests/test_integration.py::test_run_genbank --update_goldens   # one case
```

`--update_goldens` copies the freshly produced outputs in as the new reference
(applying the same skip rules) instead of asserting. Review the git diff before
committing.

## Determinism & version coupling

The committed goldens were generated on **GPU** (Pawsey, `pytest --gpu_available`),
which is also where CI validates them — both use the tolerant (drop-volatile /
≥95 % identity) path, so they are self-consistent. The same goldens also pass a
local CPU run: CPU is deterministic (`strict`), and the dropped/identity columns
mean GPU-vs-CPU differences never reach the comparison.

The goldens are coupled to:

- **phold DB version `1.0.0`** (the standard `phold install` database) — the
  Foldseek hits depend on DB content, so a DB bump requires regeneration.
- pinned dependency versions baked into the GenBank `/source` qualifier
  (e.g. `Pyrodigal-gv`) — upgrading those requires regeneration.

If a golden test starts failing only because of an intended database/dependency
update, regenerate with `--update_goldens` and commit the new references.
