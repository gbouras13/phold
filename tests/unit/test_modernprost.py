"""
Unit tests for the ModernProst (3Di + 12-state) wiring in phold.

Covers the three decision points that route a run down the ModernProst path:
``--model`` resolution, the prediction-mode auto-detection ``phold compare``
uses when no explicit flags are passed, and the Foldseek search flags. None of
these need a model, a GPU or the foldseek binary.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from phold import MODEL_CHOICES, PROSTT5_MODEL, _resolved_task, resolve_model
from phold.databases.db import check_db_12st_support
from phold.features.predict_3di_12st import mean_probs_filename
from phold.features.run_foldseek import run_foldseek_search
from phold.subcommands.compare import detect_prediction_mode


# ===========================================================================
# --model resolution
# ===========================================================================

class TestResolveModel:
    def test_prostt5_is_the_default_choice(self):
        assert MODEL_CHOICES[0] == PROSTT5_MODEL

    def test_prostt5_returns_the_cnn_checkpoint(self):
        name, model_dir, ckpt, modernprost = resolve_model(
            "prostt5", Path("/db"), finetune=False, vanilla=False
        )
        assert name == "Rostlab/ProstT5_fp16"
        assert model_dir == Path("/db")
        assert ckpt is not None and ckpt.name == "model.pt"
        assert modernprost is False

    def test_finetune_switches_encoder_and_head(self):
        name, _, ckpt, _ = resolve_model(
            "prostt5", Path("/db"), finetune=True, vanilla=False
        )
        assert name == "gbouras13/ProstT5Phold"
        assert ckpt.name == "phold_db_model.pth"

    def test_vanilla_only_switches_the_head(self):
        name, _, ckpt, _ = resolve_model(
            "prostt5", Path("/db"), finetune=True, vanilla=True
        )
        assert name == "gbouras13/ProstT5Phold"
        assert ckpt.name == "vanilla_model.pth"

    @pytest.mark.parametrize(
        "choice,expected",
        [
            ("modernprost-base", "modernprost-base"),
            ("modernprost-50M", "modernprost-50M"),
            ("modernprost-pssm", "modernprost-pssm"),
            ("modernprost-50M-pssm", "modernprost-50M-pssm"),
        ],
    )
    def test_modernprost_models_resolve(self, choice, expected):
        name, model_dir, ckpt, modernprost = resolve_model(
            choice, Path("/db"), finetune=False, vanilla=False
        )
        assert name == expected
        assert modernprost is True
        # There is no separate CNN head — the 3Di/12st heads are in the model.
        assert ckpt is None
        # The phold database directory doubles as the HuggingFace cache.
        assert model_dir == Path("/db")

    def test_lowercased_choice_maps_back_to_the_registry_spelling(self):
        # click.Choice(case_sensitive=False) lower-cases the value, but the
        # HuggingFace repo is "modernprost-50M" with a capital M.
        name, _, _, _ = resolve_model(
            "modernprost-50m", Path("/db"), finetune=False, vanilla=False
        )
        assert name == "modernprost-50M"

    def test_finetune_is_ignored_for_modernprost(self):
        name, _, ckpt, modernprost = resolve_model(
            "modernprost-base", Path("/db"), finetune=True, vanilla=True
        )
        assert name == "modernprost-base"
        assert ckpt is None
        assert modernprost is True


class TestResolvedTask:
    def test_auto_uses_the_checkpoints_trained_task(self):
        assert _resolved_task("modernprost-base", "auto") == "classification"
        assert _resolved_task("modernprost-50M", "auto") == "classification"
        assert _resolved_task("modernprost-pssm", "auto") == "pssm"
        assert _resolved_task("modernprost-50M-pssm", "auto") == "pssm"

    def test_explicit_task_overrides_the_default(self):
        assert _resolved_task("modernprost-pssm", "classification") == "classification"
        assert _resolved_task("modernprost-50M", "pssm") == "pssm"


# ===========================================================================
# Prediction-mode auto-detection
# ===========================================================================

class TestDetectPredictionMode:
    def test_prostt5_output_detects_as_3di_only(self, tmp_path):
        (tmp_path / "phold_3di.fasta").write_text(">a\nDD\n")
        assert detect_prediction_mode(tmp_path, "phold") == (False, False)

    def test_12st_fasta_detects_the_combined_alphabet(self, tmp_path):
        (tmp_path / "phold_3di.fasta").write_text(">a\nDD\n")
        (tmp_path / "phold_12st.fasta").write_text(">a\nAC\n")
        assert detect_prediction_mode(tmp_path, "phold") == (True, False)

    def test_profile_db_detects_profiles_and_implies_12st(self, tmp_path):
        profile_dir = tmp_path / "query_profiledb"
        profile_dir.mkdir()
        (profile_dir / "phold_profile_ss").write_bytes(b"\x00" * 25)
        assert detect_prediction_mode(tmp_path, "phold") == (True, True)

    def test_empty_directory_detects_as_3di_only(self, tmp_path):
        assert detect_prediction_mode(tmp_path, "phold") == (False, False)

    def test_detection_is_prefix_scoped(self, tmp_path):
        (tmp_path / "other_12st.fasta").write_text(">a\nAC\n")
        assert detect_prediction_mode(tmp_path, "phold") == (False, False)


# ===========================================================================
# Output filenames
# ===========================================================================

class TestOutputNames:
    def test_modernprost_mean_probs_filename_is_distinct(self):
        # ProstT5 and ModernProst outputs must be able to coexist in one
        # directory, and compare.py picks between them by name.
        assert mean_probs_filename("phold") == (
            "phold_modernprost_3di_mean_probabilities.csv"
        )
        assert mean_probs_filename("phold") != "phold_prostT5_3di_mean_probabilities.csv"


# ===========================================================================
# 12-state database support marker
# ===========================================================================

class TestCheckDb12stSupport:
    def test_absent_marker_reads_as_unsupported(self, tmp_path):
        assert check_db_12st_support(tmp_path) is False

    def test_marker_present(self, tmp_path):
        (tmp_path / "all_phold_structures_ss12.marker").write_text("1.1.0\n")
        assert check_db_12st_support(tmp_path) is True

    def test_directory_named_like_the_marker_does_not_count(self, tmp_path):
        (tmp_path / "all_phold_structures_ss12.marker").mkdir()
        assert check_db_12st_support(tmp_path) is False


# ===========================================================================
# Foldseek search flags
# ===========================================================================

def _search_cmd(**overrides) -> str:
    """Run run_foldseek_search with ExternalTool stubbed; return the command."""
    captured = {}

    class FakeTool:
        def __init__(self, **kwargs):
            captured["params"] = kwargs["params"]

        @staticmethod
        def run_tool(tool):
            return None

    kwargs = dict(
        query_db=Path("Q"),
        target_db=Path("T"),
        result_db=Path("R"),
        temp_db=Path("tmp"),
        threads=8,
        logdir=Path("log"),
        evalue=1e-3,
        sensitivity=9.5,
        max_seqs=1000,
        ultra_sensitive=False,
        extra_foldseek_params=None,
        foldseek_gpu=False,
        structures=False,
        clustered_db=False,
    )
    kwargs.update(overrides)

    with patch("phold.features.run_foldseek.ExternalTool", FakeTool):
        run_foldseek_search(**kwargs)

    return captured["params"]


class TestFoldseekSearchFlags:
    def test_prostt5_path_is_unchanged(self):
        cmd = _search_cmd()
        assert "--ss-12st" not in cmd
        assert "--sort-by-structure-bits" not in cmd

    def test_ss_12st_flag_added(self):
        assert "--ss-12st 1" in _search_cmd(ss_12st=True)

    def test_profiles_disable_structure_bits_sorting(self):
        cmd = _search_cmd(ss_12st=True, profiles=True)
        assert "--ss-12st 1" in cmd
        assert "--sort-by-structure-bits 0" in cmd

    def test_profiles_alone_do_not_add_ss_12st(self):
        # subcommand_compare rejects this combination before it gets here; the
        # flag builder should not paper over it by inferring one from the other.
        cmd = _search_cmd(profiles=True)
        assert "--ss-12st" not in cmd

    def test_ss_12st_survives_gpu_mode(self):
        cmd = _search_cmd(ss_12st=True, foldseek_gpu=True)
        assert "--gpu 1" in cmd
        assert "--ss-12st 1" in cmd

    def test_ss_12st_survives_ultra_sensitive(self):
        cmd = _search_cmd(ss_12st=True, ultra_sensitive=True)
        assert "--exhaustive-search" in cmd
        assert "--ss-12st 1" in cmd

    def test_extra_params_come_after_the_12st_flag(self):
        cmd = _search_cmd(ss_12st=True, extra_foldseek_params="--alignment-mode 3")
        assert cmd.index("--ss-12st 1") < cmd.index("--alignment-mode 3")
