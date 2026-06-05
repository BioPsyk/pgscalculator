#!/usr/bin/env bash
# Phase 10 acceptance test (§6.5): incremental run on REAL data, chr22.
#
#   Day 1: sBayesR only  -> prep + sumstat/weights/score/finalize (--methods sbayesr)
#   Day 2: add LDpred2   -> ldpred2 prep + weights/score/finalize (--methods ldpred2),
#                           SAME outdir, Day-1 artifacts untouched.
#
# Asserts:
#   - scores_sbayesr.gz and augmented_sbayesr.gz are byte-identical Day1 vs Day2
#   - Day 2 adds scores_ldpred2.gz + augmented_ldpred2.gz (with postp_ldpred2)
#   - per-method bench scores present for both
#
# Runs interactively in-container (no SLURM). ~10-15 min on chr22.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGS_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"
WRAPPER="${PGS_DIR}/pgscalculator-v2.sh"

DEFAULT_SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"
SUMSTAT_DIR="${SUMSTAT_DIR:-$DEFAULT_SUMSTAT_DIR}"

OUTDIR="${SCRIPT_DIR}/out_incremental_chr22"
BASE_CONFIG="${SCRIPT_DIR}/config.integration_chr22.yaml"
CONFIG_DAY1="${SCRIPT_DIR}/.config.incremental_day1.yaml"
CONFIG_DAY2="${SCRIPT_DIR}/.config.incremental_day2chr22.yaml"

# Day-1 config: same outdir, but the ldpred2: block stripped so that the prep
# group genuinely skips the LDpred2 steps (opt-in is ld_dir presence). Day 2
# uses the full config (ldpred2 enabled) pointed at the SAME outdir.
awk '
  /^ldpred2:/ {skip=1; next}
  skip && /^[^[:space:]]/ {skip=0}
  skip {next}
  {print}
' "$BASE_CONFIG" | sed "s|^outdir:.*|outdir: ${OUTDIR}|" > "$CONFIG_DAY1"
sed "s|^outdir:.*|outdir: ${OUTDIR}|" "$BASE_CONFIG" > "$CONFIG_DAY2"

echo "==================================================================="
echo "Incremental acceptance test (chr22)"
echo "  outdir : ${OUTDIR}"
echo "  sumstat: ${SUMSTAT_DIR}"
echo "==================================================================="
rm -rf "$OUTDIR"

# ---- Day 1: sBayesR only --------------------------------------------------
echo ">> Day 1: prep (sBayesR only)"
"$WRAPPER" --config "$CONFIG_DAY1" --steps prep --methods sbayesr
echo ">> Day 1: sumstat,weights,score,finalize (sBayesR only)"
"$WRAPPER" --config "$CONFIG_DAY1" --steps sumstat,weights,score,finalize --methods sbayesr -i "$SUMSTAT_DIR"

SS_OUT=$(ls -d "${OUTDIR}"/sumstats/*/ | head -1)
echo "Day-1 sumstat output: ${SS_OUT}"
[[ -f "${SS_OUT}/scores_sbayesr.gz" ]]    || { echo "FAIL: Day1 scores_sbayesr.gz missing"; exit 1; }
[[ -f "${SS_OUT}/augmented_sbayesr.gz" ]] || { echo "FAIL: Day1 augmented_sbayesr.gz missing"; exit 1; }

day1_scores_md5=$(md5sum "${SS_OUT}/scores_sbayesr.gz" | cut -d' ' -f1)
day1_aug_md5=$(md5sum "${SS_OUT}/augmented_sbayesr.gz" | cut -d' ' -f1)
echo "Day-1 scores_sbayesr.gz md5=${day1_scores_md5}"
echo "Day-1 augmented_sbayesr.gz md5=${day1_aug_md5}"

[[ -f "${SS_OUT}/augmented_ldpred2.gz" ]] && { echo "FAIL: augmented_ldpred2.gz should not exist after Day1"; exit 1; }

# ---- Day 2: add LDpred2 ----------------------------------------------------
echo ">> Day 2: prep (LDpred2 prep steps; genotypes/ldref reused)"
"$WRAPPER" --config "$CONFIG_DAY2" --steps prep --methods ldpred2
echo ">> Day 2: sumstat,weights,score,finalize (LDpred2 only)"
# Include 'sumstat' so filter-variants runs for ldpred2 (its LD-ref variant set
# differs from sBayesR, so filtered_ldpred2/ must be built). format-sumstat is
# method-agnostic and already complete from Day 1, so it is skipped.
"$WRAPPER" --config "$CONFIG_DAY2" --steps sumstat,weights,score,finalize --methods ldpred2 -i "$SUMSTAT_DIR"

# ---- Assertions ------------------------------------------------------------
echo ""
echo "==================================================================="
echo "Assertions"
echo "==================================================================="
fail=0
chk() { if eval "$2"; then echo "- [OK] $1"; else echo "- [FAIL] $1"; fail=$((fail+1)); fi; }

day2_scores_md5=$(md5sum "${SS_OUT}/scores_sbayesr.gz" | cut -d' ' -f1)
day2_aug_md5=$(md5sum "${SS_OUT}/augmented_sbayesr.gz" | cut -d' ' -f1)

chk "scores_sbayesr.gz byte-identical after Day2"    "[[ '${day1_scores_md5}' == '${day2_scores_md5}' ]]"
chk "augmented_sbayesr.gz byte-identical after Day2" "[[ '${day1_aug_md5}' == '${day2_aug_md5}' ]]"
chk "scores_ldpred2.gz created"                      "[[ -f '${SS_OUT}/scores_ldpred2.gz' ]]"
chk "augmented_ldpred2.gz created"                   "[[ -f '${SS_OUT}/augmented_ldpred2.gz' ]]"
chk "bench_score_sbayesr.gz present"                 "[[ -f '${SS_OUT}/bench_score_sbayesr.gz' ]]"
chk "bench_score_ldpred2.gz present"                 "[[ -f '${SS_OUT}/bench_score_ldpred2.gz' ]]"
chk "ldpred2 diagnostics summary.tsv present"        "[[ -f '${SS_OUT}/details/ldpred2/summary.tsv' ]]"

if [[ -f "${SS_OUT}/augmented_ldpred2.gz" ]]; then
  # sed -n '1p' consumes the whole stream (no early close), so zcat is not
  # killed by SIGPIPE -- which under `set -o pipefail` would otherwise abort.
  hdr=$(zcat "${SS_OUT}/augmented_ldpred2.gz" | sed -n '1p')
  chk "augmented_ldpred2 header has postEffect"  "[[ '$hdr' == *postEffect* ]]"
  chk "augmented_ldpred2 header has postp_ldpred2" "[[ '$hdr' == *postp_ldpred2* ]]"
  chk "augmented_ldpred2 header has benchEffect" "[[ '$hdr' == *benchEffect* ]]"
fi

echo ""
if [[ $fail -gt 0 ]]; then
  echo "INCREMENTAL ACCEPTANCE: FAILED (${fail})"
  exit 1
fi
echo "INCREMENTAL ACCEPTANCE: PASSED"
