#!/usr/bin/env bash
# run_comparison.sh
# Runs all output-producing phold tests with both:
#   - pholdENV  (dev: phold_lib branch with pholdlib)
#   - phold_125 (bioconda v1.2.5 reference)
# then compares key output files ignoring timestamps.
#
# Usage (from repo root):
#   bash tests/run_comparison.sh
#
# Outputs go to /tmp/phold_compare/

set -uo pipefail

# Log any unexpected exit with the line number so we can diagnose why the
# script dies mid-run.
trap 'echo "[EXIT] script exiting: line=$LINENO status=$? signal=${_last_sig:-none}" >> /tmp/phold_compare_exit.log' EXIT
trap '_last_sig=URG'  URG
trap '_last_sig=HUP'  HUP
trap '_last_sig=TERM' TERM
trap '_last_sig=INT'  INT

# Resolve conda binary in this order: $CONDA_EXE → conda on PATH →
# common miniforge/miniconda install locations. Daemonised shells (no
# login profile) often have none of the first two set.
_resolve_conda() {
    if [ -n "${CONDA_EXE:-}" ] && [ -x "$CONDA_EXE" ]; then echo "$CONDA_EXE"; return; fi
    local found="$(command -v conda 2>/dev/null)"
    if [ -n "$found" ]; then echo "$found"; return; fi
    for candidate in \
        "$HOME/miniforge3/bin/conda" \
        "$HOME/miniconda3/bin/conda" \
        "$HOME/anaconda3/bin/conda" \
        "/opt/miniforge3/bin/conda" \
        "/opt/conda/bin/conda"; do
        if [ -x "$candidate" ]; then echo "$candidate"; return; fi
    done
    echo "" # nothing found — error below
}
CONDA="$(_resolve_conda)"
if [ -z "$CONDA" ]; then
    echo "ERROR: could not find a conda binary. Set CONDA_EXE or add conda to PATH." >&2
    exit 2
fi
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="${REPO}/tests/test_data/phold_db"
TD="${REPO}/tests/test_data"
BASE="/tmp/phold_compare"
T=1   # single thread for determinism

DEV_ENV="pholdENV"
REF_ENV="phold_125"

DEV_OUT="${BASE}/dev_outputs"
REF_OUT="${BASE}/ref_outputs"

mkdir -p "$BASE" "$DEV_OUT" "$REF_OUT"

# CSV of per-case wall-clock per environment; appended to as we go.
TIMING_CSV="${BASE}/timings.csv"
echo "case,env,wall_seconds,status" > "$TIMING_CSV"

PASS=0
FAIL=0
ERRORS=()

# ── helpers ──────────────────────────────────────────────────────────────────

run_test() {
    local env="$1"
    local label="$2"
    local log="$3"
    shift 3
    local cmd="$*"

    # Bypass ``conda run`` and call the env's binaries directly. ``conda
    # run`` is heavyweight (spawns a Python interpreter just to set up
    # env vars) and was observed to die silently when daemonised.
    # Resolve <conda_root>/envs/<env>/bin and prepend it to PATH for
    # the child shell — that's exactly what ``conda activate`` does.
    local conda_root="$(dirname "$(dirname "$CONDA")")"
    local env_bin="${conda_root}/envs/${env}/bin"
    if [ ! -d "$env_bin" ]; then
        echo "ERROR: env bin dir not found: $env_bin" >&2
        return 2
    fi

    printf "  [%s] %s starting...\n" "$(date '+%H:%M:%S')" "$label"
    local t0=$SECONDS
    local rc=0
    if PATH="${env_bin}:${PATH}" bash -c "cd '$REPO' && $cmd" > "$log" 2>&1; then
        local dt=$((SECONDS - t0))
        printf "  [%s] %s DONE ✓ (%ds)\n" "$(date '+%H:%M:%S')" "$label" "$dt"
        # Extract case name and env from label "case[env]" → "case,env"
        local case_name="${label%\[*}"
        local env_tag="${label#*\[}"; env_tag="${env_tag%\]}"
        echo "$case_name,$env_tag,$dt,ok" >> "$TIMING_CSV"
        return 0
    else
        rc=$?
        local dt=$((SECONDS - t0))
        printf "  [%s] %s FAILED ✗ (%ds, see %s)\n" "$(date '+%H:%M:%S')" "$label" "$dt" "$log"
        local case_name="${label%\[*}"
        local env_tag="${label#*\[}"; env_tag="${env_tag%\]}"
        echo "$case_name,$env_tag,$dt,fail" >> "$TIMING_CSV"
        ERRORS+=("$label FAILED — see $log")
        FAIL=$((FAIL + 1))
        return 1
    fi
}

run_case() {
    local name="$1"; shift
    local phold_cmd="$*"

    printf "\n══ %-50s ══\n" "$name"

    local dev_log="${BASE}/${name}_dev.log"
    local ref_log="${BASE}/${name}_ref.log"
    local out_dir="${TD}/outputs"

    # ── dev run ────────────────────────────────────────────
    rm -rf "$out_dir"
    mkdir -p "$out_dir"

    if run_test "$DEV_ENV" "${name}[dev]" "$dev_log" "$phold_cmd"; then
        # snapshot dev outputs
        local snap_dev="${DEV_OUT}/${name}"
        rm -rf "$snap_dev"
        cp -r "$out_dir" "$snap_dev"
    else
        return
    fi

    # ── ref run ────────────────────────────────────────────
    rm -rf "$out_dir"
    mkdir -p "$out_dir"

    if run_test "$REF_ENV" "${name}[ref]" "$ref_log" "$phold_cmd"; then
        local snap_ref="${REF_OUT}/${name}"
        rm -rf "$snap_ref"
        cp -r "$out_dir" "$snap_ref"
    else
        return
    fi

    # ── compare ────────────────────────────────────────────
    printf "  Comparing outputs...\n"
    if python3 "${REPO}/tests/compare_outputs.py" "${DEV_OUT}/${name}" "${REF_OUT}/${name}" --cpu; then
        printf "  ✓ MATCH\n"
        PASS=$((PASS + 1))
    else
        printf "  ✗ DIFFER\n"
        ERRORS+=("$name: outputs differ")
        FAIL=$((FAIL + 1))
    fi
}

# ── test cases ───────────────────────────────────────────────────────────────
# predict: pure 3Di inference (most critical — tests pholdlib refactor)

run_case "predict_genbank" \
    "phold predict -i ${TD}/combined_truncated_acr_defense_vfdb_card.gbk \
     -o ${TD}/outputs/predict_gbk -t $T -d $DB -f --cpu"

run_case "predict_fasta" \
    "phold predict -i ${TD}/combined_truncated_acr_defense_vfdb_card.fasta \
     -o ${TD}/outputs/predict_fasta -t $T -d $DB -f --cpu"

run_case "proteins_predict" \
    "phold proteins-predict -i ${TD}/phanotate.faa \
     -o ${TD}/outputs/proteins_predict -t $T -d $DB -f --cpu"

run_case "proteins_predict_gzip" \
    "phold proteins-predict -i ${TD}/phanotate.faa.gz \
     -o ${TD}/outputs/proteins_predict_gz -t $T -d $DB -f --cpu"

# run: full pipeline (predict + foldseek)

run_case "run_genbank" \
    "phold run -i ${TD}/NC_043029_pharokka1.4.1.gbk \
     -o ${TD}/outputs/run_gbk -t $T -d $DB -f --cpu"

run_case "run_fasta" \
    "phold run -i ${TD}/combined_truncated_acr_defense_vfdb_card.fasta \
     -o ${TD}/outputs/run_fasta -t $T -d $DB -f --cpu"

run_case "run_genbank_ncbi" \
    "phold run -i ${TD}/NC_043029_ncbi.gbk \
     -o ${TD}/outputs/run_gbk_ncbi -t $T -d $DB -f --cpu"

run_case "run_genbank_bakta" \
    "phold run -i ${TD}/NC_043029_bakta.gbk \
     -o ${TD}/outputs/run_gbk_bakta -t $T -d $DB -f --cpu"

run_case "run_genbank_old_pharokka" \
    "phold run -i ${TD}/combined_truncated_acr_defense_vfdb_card.gbk \
     -o ${TD}/outputs/run_gbk_old -t $T -d $DB -f --cpu"

run_case "run_genbank_long_header" \
    "phold run -i ${TD}/long_header.gbk \
     -o ${TD}/outputs/run_gbk_long -t $T -d $DB -f --cpu"

run_case "run_fasta_long_header" \
    "phold run -i ${TD}/long_header.fasta \
     -o ${TD}/outputs/run_fasta_long -t $T -d $DB -f --cpu"

run_case "run_efam" \
    "phold run -i ${TD}/KF623293.1_subset_efam.fasta \
     -o ${TD}/outputs/run_efam -t $T -d $DB -f --cpu"

run_case "run_netflax" \
    "phold run -i ${TD}/WP_006719989_subset_test.fasta \
     -o ${TD}/outputs/run_netflax -t $T -d $DB -f --cpu"

# compare: uses pre-built predict output

run_case "predict_then_compare_genbank" \
    "phold predict -i ${TD}/combined_truncated_acr_defense_vfdb_card.gbk \
       -o ${TD}/outputs/pred_for_cmp -t $T -d $DB -f --cpu && \
     phold compare -i ${TD}/combined_truncated_acr_defense_vfdb_card.gbk \
       --predictions_dir ${TD}/outputs/pred_for_cmp \
       -o ${TD}/outputs/compare_gbk -t $T -d $DB -f"

run_case "predict_then_compare_fasta" \
    "phold predict -i ${TD}/combined_truncated_acr_defense_vfdb_card.fasta \
       -o ${TD}/outputs/pred_fasta_for_cmp -t $T -d $DB -f --cpu && \
     phold compare -i ${TD}/combined_truncated_acr_defense_vfdb_card.fasta \
       --predictions_dir ${TD}/outputs/pred_fasta_for_cmp \
       -o ${TD}/outputs/compare_fasta -t $T -d $DB -f"

run_case "proteins_predict_then_compare" \
    "phold proteins-predict -i ${TD}/phanotate.faa \
       -o ${TD}/outputs/prot_pred_for_cmp -t $T -d $DB -f --cpu && \
     phold proteins-compare -i ${TD}/phanotate.faa \
       --predictions_dir ${TD}/outputs/prot_pred_for_cmp \
       -o ${TD}/outputs/proteins_compare -t $T -d $DB -f"

# ── summary ──────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════════"
printf "  Results: %d passed, %d failed\n" "$PASS" "$FAIL"
echo "══════════════════════════════════════════════════════"
if [ "${#ERRORS[@]}" -gt 0 ]; then
    echo "  Failures:"
    for e in "${ERRORS[@]}"; do
        echo "    ✗ $e"
    done
    exit 1
else
    echo "  All cases match ✓"
    exit 0
fi
