#!/usr/bin/env bash
# Unit tests for Phase 8 driver helpers (defined in pgscalculator-v2.sh).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_DIR}/bin/lib/common.sh"

# Re-define driver helpers (kept in sync with pgscalculator-v2.sh).
build_sumstat_filter_chr_cmds() {
  local run_cmd="$1"
  local body=""
  local m
  for m in ${CFG_METHODS:-sbayesr}; do
    body="${body}${run_cmd} --method ${m} --_chr \"\$CHR\"; "
  done
  echo "$body"
}

score_profile_mapped_dir() {
  local outdir_host="$1"
  local sumstat_name="$2"
  local step_profile="$3"
  local method=""
  case "$step_profile" in
    score_sbayesr|weights_sbayesr) method="sbayesr" ;;
    score_ldpred2) method="ldpred2" ;;
    score) method="sbayesr" ;;
    *) echo ""; return 0 ;;
  esac
  local base="${outdir_host}/sumstats/${sumstat_name}/work"
  if [[ "$method" == "sbayesr" ]]; then
    if [[ -d "${base}/posteriors_mapped" ]]; then
      echo "${base}/posteriors_mapped"
    fi
  else
    echo "${base}/posteriors_mapped_ldpred2"
  fi
}

VERBOSE=0
failures=0
tmpdir=""

cleanup() {
  [[ -n "$tmpdir" && -d "$tmpdir" ]] && rm -rf "$tmpdir"
}
trap cleanup EXIT

assert_contains() {
  local label="$1"
  local haystack="$2"
  local needle="$3"
  if [[ "$haystack" != *"$needle"* ]]; then
    echo "- [FAIL] ${label}: expected to contain '${needle}'"
    failures=$((failures + 1))
  else
    echo "- [OK] ${label}"
  fi
}

echo ">> Test driver dispatch helpers"

tmpdir=$(mktemp -d)
export CFG_METHODS="sbayesr ldpred2"

body=$(build_sumstat_filter_chr_cmds "RUN")
assert_contains "filter sbayesr" "$body" "--method sbayesr"
assert_contains "filter ldpred2" "$body" "--method ldpred2"

mkdir -p "${tmpdir}/sumstats/trait1/work/posteriors_mapped_ldpred2"
dir=$(score_profile_mapped_dir "$tmpdir" "trait1" "score_ldpred2")
if [[ "$dir" == "${tmpdir}/sumstats/trait1/work/posteriors_mapped_ldpred2" ]]; then
  echo "- [OK] score_ldpred2 mapped dir"
else
  echo "- [FAIL] score_ldpred2 mapped dir: got '${dir}'"
  failures=$((failures + 1))
fi

if [[ $failures -gt 0 ]]; then
  echo "Tests failed: ${failures}"
  exit 1
fi
echo "All driver dispatch tests passed."
