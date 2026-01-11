#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGS_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"

WRAPPER="${PGS_DIR}/pgscalculator-v2.sh"
CONFIG="${SCRIPT_DIR}/config.yaml"

DEFAULT_SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"
SUMSTAT_DIR="${SUMSTAT_DIR:-$DEFAULT_SUMSTAT_DIR}"
if [[ -z "$SUMSTAT_DIR" ]]; then
  >&2 echo "Error: set SUMSTAT_DIR to a cleansumstats output folder, e.g.:"
  >&2 echo "  export SUMSTAT_DIR=/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"
  exit 1
fi
if [[ ! -d "$SUMSTAT_DIR" ]]; then
  >&2 echo "Error: SUMSTAT_DIR does not exist: ${SUMSTAT_DIR}"
  exit 1
fi

echo "Config: ${CONFIG}"
echo "Sumstat: ${SUMSTAT_DIR}"
echo ""

echo "== Running: prep (interactive) =="
"$WRAPPER" --config "$CONFIG" --steps prep

echo ""
echo "== Running: sumstat,posteriors,score (interactive) =="
"$WRAPPER" --config "$CONFIG" --steps sumstat,posteriors,score -i "$SUMSTAT_DIR"

echo ""
echo "Done. Outdir:"
grep -E '^outdir:' "$CONFIG" || true

