#!/usr/bin/env bash
# Unit tests for Phase 9 discovery helpers and combine-scores naming.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_DIR}/bin/lib/common.sh"

VERBOSE=0
tmpdir=""
failures=0

cleanup() {
    [[ -n "$tmpdir" && -d "$tmpdir" ]] && rm -rf "$tmpdir"
}
trap cleanup EXIT

assert_eq() {
    local label="$1" got="$2" want="$3"
    if [[ "$got" != "$want" ]]; then
        echo "- [FAIL] ${label}: got '${got}', want '${want}'"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

echo ">> Test Phase 9 discovery helpers"
tmpdir=$(mktemp -d)
sumstat_dir="${tmpdir}/sumstats/trait1"
work="${sumstat_dir}/work"
mkdir -p "${work}/posteriors_mapped" "${work}/scores" "${work}/filtered_sbayesr"

echo -e "ID\tA1\tA2\tFreq\tEffect\tSE\tPIP" > "${work}/posteriors_mapped/chr22.snpRes"
echo "rs1 A G 0.5 0.01 NA 0.5" >> "${work}/posteriors_mapped/chr22.snpRes"
echo -e "FID\tIID\tSCORE1_SUM" > "${work}/scores/chr22.sscore"
echo "0 sample1 1.0" >> "${work}/scores/chr22.sscore"
echo -e "CHR\tPOS" > "${work}/filtered_sbayesr/chr22_matched.tsv"
echo "22 1" >> "${work}/filtered_sbayesr/chr22_matched.tsv"

discovered=$(discover_posterior_methods "$sumstat_dir")
assert_eq "discover posteriors sbayesr" "$discovered" "sbayesr"

discovered=$(discover_score_methods "$sumstat_dir")
assert_eq "discover scores sbayesr" "$discovered" "sbayesr"

assert_eq "scores gz name" "$(method_scores_gz_name sbayesr)" "scores_sbayesr.gz"
assert_eq "ldpred2 gz name" "$(method_scores_gz_name ldpred2)" "scores_ldpred2.gz"

mkdir -p "${work}/posteriors_mapped_ldpred2"
echo -e "ID\tA1\tA2\tFreq\tEffect\tSE\tPIP" > "${work}/posteriors_mapped_ldpred2/chr22.snpRes"
echo "rs1 A G 0.5 0.02 NA 0.6" >> "${work}/posteriors_mapped_ldpred2/chr22.snpRes"
export CFG_METHODS="sbayesr ldpred2"
ordered=$(order_discovered_methods "$(discover_posterior_methods "$sumstat_dir")")
assert_eq "ordered both methods" "$ordered" "sbayesr ldpred2"

if [[ $failures -gt 0 ]]; then
    echo "Tests failed: ${failures}"
    exit 1
fi
echo "All Phase 9 discovery tests passed."
