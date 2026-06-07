#!/usr/bin/env bash
# Unit test: format-posteriors LDpred2 path (12-col variant_map_ldpred2 + LDpred2 .snpRes layout)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
export SCRIPT_DIR="${PROJECT_DIR}/bin"
export LIB_DIR="${PROJECT_DIR}/bin/lib"
export STEPS_DIR="${LIB_DIR}/steps"
VERBOSE=0
DRY_RUN=0

source "${LIB_DIR}/common.sh"
source "${STEPS_DIR}/format_posteriors.sh"

tmpdir=""
failures=0

cleanup() {
    [[ -n "$tmpdir" && -d "$tmpdir" ]] && rm -rf "$tmpdir"
}
trap cleanup EXIT

assert_gt_zero() {
    local label="$1"
    local val="$2"
    if [[ ! "$val" =~ ^[0-9]+$ ]] || [[ "$val" -le 0 ]]; then
        echo "- [FAIL] ${label}: expected mapped count > 0, got '${val}'"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

echo ">> Test format-posteriors LDpred2 mapping"

tmpdir=$(mktemp -d)

# 12-col LDpred2 variant map (cols 4=geno_snpid, 7=ldref_snpid)
cat > "${tmpdir}/variant_map_ldpred2.tsv" <<'EOF'
chr	pos_b37	pos_b38	geno_snpid	geno_a1	geno_a2	ldref_snpid	ldref_a1	ldref_a2	ldref_a2freq	ld	block_id
22	17309296	16828406	22:16828406_A_G	A	G	rs175146	A	G	0.11	0.5	1
22	18181984	17699218	22:17699218_T_G	T	G	rs17207051	T	G	0.02	0.4	2
EOF

# LDpred2 calc-ldpred2 .snpRes layout (Name = ldref_snpid in column 2)
cat > "${tmpdir}/chr22.snpRes" <<'EOF'
Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP LastSampleEff
1 rs175146 22 17309296 A G 0.890000 -0.001149 NA 0.61250019 -0.001149
2 rs17207051 22 18181984 T G 0.980000 -0.001620 NA 0.65374994 -0.001620
EOF

rsid_map="${tmpdir}/ldref_to_genoid.tsv"
ensure_rsid_mapping "${tmpdir}/variant_map_ldpred2.tsv" "$rsid_map"

out="${tmpdir}/chr22_mapped.snpRes"
mapped=$(map_posteriors_for_chr 22 "${tmpdir}/chr22.snpRes" "$rsid_map" "$out")

assert_gt_zero "mapped variant count" "$mapped"

if ! awk -F'\t' '$1=="22:16828406_A_G"{found=1} END{exit !found}' "$out"; then
    echo "- [FAIL] output missing 22:16828406_A_G genotype row"
    failures=$((failures + 1))
else
    echo "- [OK] 22:16828406_A_G present in mapped output"
fi

if [[ $failures -gt 0 ]]; then
    echo "Tests failed: ${failures}"
    exit 1
fi
echo "All format-posteriors LDpred2 tests passed."
