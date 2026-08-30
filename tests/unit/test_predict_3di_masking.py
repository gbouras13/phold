"""
Regression tests for 3Di masking when per-residue probabilities are absent.

``--omit_probs`` leaves ``all_prob`` as ``None`` in the
``(pred, mean_prob, all_prob)`` tuple. ``write_predictions`` used to index
``all_prob[0]`` unconditionally — including when masking was disabled
(``mask_threshold=0``) — so the run died with
``TypeError: 'NoneType' object is not subscriptable`` at write time, after
the entire GPU prediction had already completed.

Run::

    pytest tests/unit/test_predict_3di_masking.py -v
"""
from __future__ import annotations

import numpy as np
import pytest

from phold.features.predict_3Di import write_predictions
from phold.utils.validation import validate_mask_options


def _pred(labels):
    """A prediction array in the np.byte form the CNN head emits."""
    return np.array(labels, dtype=np.byte)


def _read_fasta(path):
    """Tiny FASTA reader -> {header: sequence}."""
    records, header = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:]
            records[header] = ""
        elif header is not None:
            records[header] += line
    return records


# ─── write_predictions: masking disabled ────────────────────────────────────


def test_write_predictions_mask_threshold_zero_with_none_probs(tmp_path):
    """mask_threshold=0 + all_prob=None — the original crash.

    Masking is off, so the probabilities are never needed and the loop must
    not be entered at all.
    """
    predictions = {"contig1": {"cds1": (_pred([0, 1, 2]), 82.69, None)}}
    out = tmp_path / "out_3di.fasta"

    write_predictions(predictions, out, proteins_flag=False, mask_threshold=0)

    assert _read_fasta(out) == {"contig1:cds1": "ACD"}


def test_write_predictions_mask_threshold_zero_leaves_prediction_untouched(tmp_path):
    """mask_threshold=0 must not mask, even when probabilities are available.

    ``all_prob[0] < 0`` is all-False, so the old loop was a no-op that still
    required the array; the output is unchanged but nothing is now indexed.
    """
    pred = _pred([0, 1, 2])
    predictions = {"contig1": {"cds1": (pred, 82.69, np.array([[0.0, 0.0, 0.0]], dtype=np.float32))}}
    out = tmp_path / "out_3di.fasta"

    write_predictions(predictions, out, proteins_flag=False, mask_threshold=0)

    assert _read_fasta(out) == {"contig1:cds1": "ACD"}
    assert pred.tolist() == [0, 1, 2]


# ─── write_predictions: masking enabled ─────────────────────────────────────


def test_write_predictions_masking_skipped_when_probs_missing(tmp_path):
    """mask_threshold>0 + all_prob=None — skip with a warning, never crash.

    The CLI rejects this combination up front (see
    ``test_validate_mask_options_*``); this covers library callers of
    ``get_embeddings``/``write_predictions`` directly.
    """
    predictions = {"contig1": {"cds1": (_pred([0, 1, 2]), 82.69, None)}}
    out = tmp_path / "out_3di.fasta"

    write_predictions(predictions, out, proteins_flag=False, mask_threshold=25)

    assert _read_fasta(out) == {"contig1:cds1": "ACD"}


def test_write_predictions_masking_applied(tmp_path):
    """The masking path itself still works: probs below the threshold -> 'X'."""
    all_prob = np.array([[0.9, 0.1, 0.5]], dtype=np.float32)
    predictions = {"contig1": {"cds1": (_pred([0, 1, 2]), 50.0, all_prob)}}
    out = tmp_path / "out_3di.fasta"

    # threshold 25 -> mask_prop 0.25, so only the 0.1 residue is masked
    write_predictions(predictions, out, proteins_flag=False, mask_threshold=25)

    assert _read_fasta(out) == {"contig1:cds1": "AXD"}


def test_write_predictions_mixed_none_and_present_probs(tmp_path):
    """One sequence with probs, one without — the second must not poison the first."""
    predictions = {
        "contig1": {
            "with_probs": (_pred([0, 1, 2]), 50.0, np.array([[0.9, 0.1, 0.9]], dtype=np.float32)),
            "no_probs": (_pred([0, 1, 2]), 50.0, None),
        }
    }
    out = tmp_path / "out_3di.fasta"

    write_predictions(predictions, out, proteins_flag=False, mask_threshold=25)

    assert _read_fasta(out) == {
        "contig1:with_probs": "AXD",
        "contig1:no_probs": "ACD",
    }


def test_write_predictions_proteins_flag_headers(tmp_path):
    """proteins mode: bare seq_id headers, and still no probs required."""
    predictions = {"proteins": {"prot1": (_pred([0, 1]), 82.69, None)}}
    out = tmp_path / "out_3di.fasta"

    write_predictions(predictions, out, proteins_flag=True, mask_threshold=0)

    assert _read_fasta(out) == {"prot1": "AC"}


def test_write_predictions_drops_zero_length(tmp_path):
    """Zero-length predictions are still dropped (issue #47) with all_prob=None."""
    predictions = {
        "contig1": {
            "empty": (_pred([]), 0.0, None),
            "ok": (_pred([0]), 82.69, None),
        }
    }
    out = tmp_path / "out_3di.fasta"

    write_predictions(predictions, out, proteins_flag=False, mask_threshold=0)

    assert _read_fasta(out) == {"contig1:ok": "A"}


# ─── validate_mask_options: fail fast at argument-parsing time ──────────────


@pytest.mark.parametrize("mask_threshold", [25, 0.1, 100])
def test_validate_mask_options_rejects_incompatible_combo(mask_threshold):
    """--omit_probs with masking on must exit before any model is loaded."""
    with pytest.raises(SystemExit) as exc:
        validate_mask_options(omit_probs=True, mask_threshold=mask_threshold)
    assert exc.value.code == 1


@pytest.mark.parametrize(
    "omit_probs, mask_threshold",
    [
        (True, 0),    # masking explicitly disabled — the documented workaround
        (False, 25),  # the default: probs retained, masking on
        (False, 0),
    ],
)
def test_validate_mask_options_accepts_valid_combos(omit_probs, mask_threshold):
    validate_mask_options(omit_probs=omit_probs, mask_threshold=mask_threshold)
