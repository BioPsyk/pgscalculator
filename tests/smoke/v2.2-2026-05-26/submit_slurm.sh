#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGS_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"
WRAPPER="${PGS_DIR}/pgscalculator-v2.sh"

CONFIG_NAME="${1:-config.both.yaml}"
CONFIG_SRC="${SCRIPT_DIR}/${CONFIG_NAME}"
DEFAULT_SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"
SUMSTAT_DIR="${SUMSTAT_DIR:-$DEFAULT_SUMSTAT_DIR}"
DEFAULT_LDPRED2_LD_DIR="${PGS_DIR}/references/ld-ldpred2/hm3"
LDPRED2_LD_DIR="${LDPRED2_LD_DIR:-$DEFAULT_LDPRED2_LD_DIR}"

CONFIG="${SCRIPT_DIR}/.config.resolved.yaml"
sed "s|__LDPRED2_LD_DIR__|${LDPRED2_LD_DIR}|g" "$CONFIG_SRC" > "$CONFIG"

if [[ ! -d "$SUMSTAT_DIR" ]]; then
  >&2 echo "Error: SUMSTAT_DIR does not exist: ${SUMSTAT_DIR}"
  exit 1
fi

echo "Config: ${CONFIG}"
echo "Sumstat: ${SUMSTAT_DIR}"
echo ""

echo "== Submitting: prep (SLURM) =="
"$WRAPPER" --config "$CONFIG" --steps prep --sbatch

echo ""
echo "== Submitting: sumstat,weights,score,finalize (SLURM) =="
"$WRAPPER" --config "$CONFIG" --steps sumstat,weights,score,finalize --sbatch -i "$SUMSTAT_DIR"
