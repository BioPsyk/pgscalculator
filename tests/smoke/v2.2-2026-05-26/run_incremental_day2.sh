#!/usr/bin/env bash
# Day 2 of §6.5: add LDpred2 to an existing sBayesR-only outdir (config.sbayesr.yaml outdir).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGS_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"
WRAPPER="${PGS_DIR}/pgscalculator.sh"

DEFAULT_SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"
SUMSTAT_DIR="${SUMSTAT_DIR:-$DEFAULT_SUMSTAT_DIR}"
DEFAULT_LDPRED2_LD_DIR="${PGS_DIR}/references/ld-ldpred2/hm3"
LDPRED2_LD_DIR="${LDPRED2_LD_DIR:-$DEFAULT_LDPRED2_LD_DIR}"

if [[ ! -d "$LDPRED2_LD_DIR" ]]; then
  >&2 echo "Error: set LDPRED2_LD_DIR to HM3 LD reference directory"
  exit 1
fi

# Reuse Day-1 outdir but run only LDpred2 weights/score/finalize
CONFIG_SRC="${SCRIPT_DIR}/config.ldpred2.yaml"
CONFIG="${SCRIPT_DIR}/.config.incremental_day2.yaml"
sed "s|__LDPRED2_LD_DIR__|${LDPRED2_LD_DIR}|g" "$CONFIG_SRC" > "$CONFIG"
# Point at sBayesR smoke outdir so Day-1 artifacts are preserved
sed -i 's|out_ldpred2|out_sbayesr|' "$CONFIG"

echo "Incremental Day 2 — LDpred2 only, outdir from config.sbayesr.yaml"
echo "Config: ${CONFIG}"
echo "Sumstat: ${SUMSTAT_DIR}"
echo ""

# Include 'sumstat' so filter-variants runs for ldpred2 (its LD-ref variant set
# differs from sBayesR). format-sumstat is method-agnostic and already complete
# from Day 1, so it is skipped; only filter-variants (ldpred2) is new.
"$WRAPPER" --config "$CONFIG" --steps sumstat,weights,score,finalize \
  --methods ldpred2 -i "$SUMSTAT_DIR"

echo "Done. Compare scores_sbayesr.gz to Day-1 copy; augmented_sumstat should list both postEffect_* columns."
