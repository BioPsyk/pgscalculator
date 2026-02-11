#!/usr/bin/env bash
# Smoke test: fresh outdir, all chromosomes, all steps via --sbatch.
# Submits prep, then the per-sumstat driver (sumstat + weights + score + finalize) with
# --dependency=afterok:PREP_JOBID so it runs only after prep completes.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGS_DIR="$(dirname "$(dirname "$(dirname "$SCRIPT_DIR")")")"

WRAPPER="${PGS_DIR}/pgscalculator-v2.sh"
CONFIG="${SCRIPT_DIR}/config.yaml"

# Output directory for this run (removed at start for a from-scratch run)
OUTDIR="${SCRIPT_DIR}/out_sbatch"

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

# From-scratch: remove existing outdir so prep and sumstat run fully
if [[ -d "$OUTDIR" ]]; then
  echo "Removing existing outdir for from-scratch run: ${OUTDIR}"
  rm -rf "$OUTDIR"
fi

echo "Config: ${CONFIG}"
echo "Outdir: ${OUTDIR}"
echo "Sumstat: ${SUMSTAT_DIR}"
echo ""

echo "== Submitting: prep (SLURM) =="
PREP_OUT=$("$WRAPPER" --config "$CONFIG" --steps prep --sbatch -o "$OUTDIR" 2>&1)
echo "$PREP_OUT"

PREP_JOBID=$(echo "$PREP_OUT" | sed -n 's/^Submitted: *//p' | tr -d '\r\n')
if [[ -z "$PREP_JOBID" ]]; then
  >&2 echo "Error: could not get prep job id from submission output"
  exit 1
fi
echo "Prep job id: ${PREP_JOBID}"
echo ""

echo "== Submitting: sumstat,weights,score,finalize (SLURM driver job, after prep) =="
export SLURM_DEPENDENCY="afterok:${PREP_JOBID}"
"$WRAPPER" --config "$CONFIG" --steps sumstat,weights,score,finalize --sbatch -i "$SUMSTAT_DIR" -o "$OUTDIR"
unset SLURM_DEPENDENCY

echo ""
echo "Done. Outputs will be under: ${OUTDIR}"
