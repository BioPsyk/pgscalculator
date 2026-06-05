#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGS_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"
WRAPPER="${PGS_DIR}/pgscalculator.sh"

CONFIG_NAME="${1:-config.both.yaml}"
CONFIG_SRC="${SCRIPT_DIR}/${CONFIG_NAME}"
if [[ ! -f "$CONFIG_SRC" ]]; then
  >&2 echo "Usage: $0 [config.sbayesr.yaml|config.ldpred2.yaml|config.both.yaml]"
  exit 1
fi

DEFAULT_SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"
SUMSTAT_DIR="${SUMSTAT_DIR:-$DEFAULT_SUMSTAT_DIR}"
DEFAULT_LDPRED2_LD_DIR="${PGS_DIR}/references/ld-ldpred2/hm3"
LDPRED2_LD_DIR="${LDPRED2_LD_DIR:-$DEFAULT_LDPRED2_LD_DIR}"

if [[ ! -d "$SUMSTAT_DIR" ]]; then
  >&2 echo "Error: SUMSTAT_DIR does not exist: ${SUMSTAT_DIR}"
  exit 1
fi

CONFIG="${SCRIPT_DIR}/.config.resolved.yaml"
sed "s|__LDPRED2_LD_DIR__|${LDPRED2_LD_DIR}|g" "$CONFIG_SRC" > "$CONFIG"

if grep -q '__LDPRED2_LD_DIR__' "$CONFIG_SRC" 2>/dev/null || grep -q 'ldpred2:' "$CONFIG_SRC"; then
  if [[ ! -d "$LDPRED2_LD_DIR" ]]; then
    >&2 echo "Error: LDPRED2_LD_DIR does not exist: ${LDPRED2_LD_DIR}"
    >&2 echo "Download LD reference per docs/references.md or set LDPRED2_LD_DIR."
    exit 1
  fi
fi

echo "Config: ${CONFIG} (from ${CONFIG_NAME})"
echo "Sumstat: ${SUMSTAT_DIR}"
[[ -d "$LDPRED2_LD_DIR" ]] && echo "LDpred2 LD: ${LDPRED2_LD_DIR}"
echo ""

echo "== Running: prep (interactive) =="
"$WRAPPER" --config "$CONFIG" --steps prep

echo ""
echo "== Running: sumstat,weights,score,finalize (interactive) =="
"$WRAPPER" --config "$CONFIG" --steps sumstat,weights,score,finalize -i "$SUMSTAT_DIR"

echo ""
echo "Done. Outdir:"
grep -E '^outdir:' "$CONFIG" || true
