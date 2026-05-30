#!/usr/bin/env bash
# Acceptance test for §6.5/§11 incremental runs: a Day-1 sBayesR finalize produces
# augmented_sbayesr.gz; a Day-2 LDpred2 finalize adds augmented_ldpred2.gz without
# rewriting scores_sbayesr.gz or augmented_sbayesr.gz. Each method gets its own
# self-contained augmented file (with its own postEffect/benchEffect inside).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
export SCRIPT_DIR="${PROJECT_DIR}/bin"
export LIB_DIR="${PROJECT_DIR}/bin/lib"
export STEPS_DIR="${LIB_DIR}/steps"
VERBOSE=0
DRY_RUN=0

source "${LIB_DIR}/common.sh"
source "${STEPS_DIR}/combine_scores.sh"
source "${STEPS_DIR}/finalize_output.sh"

tmpdir=""
failures=0

cleanup() {
    [[ -n "$tmpdir" && -d "$tmpdir" ]] && rm -rf "$tmpdir"
}
trap cleanup EXIT

assert_contains() {
    local label="$1" haystack="$2" needle="$3"
    if [[ "$haystack" != *"$needle"* ]]; then
        echo "- [FAIL] ${label}: expected header to contain '${needle}'"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

assert_not_contains() {
    local label="$1" haystack="$2" needle="$3"
    if [[ "$haystack" == *"$needle"* ]]; then
        echo "- [FAIL] ${label}: header should not contain '${needle}'"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

assert_eq_file() {
    local label="$1" a="$2" b="$3"
    if ! cmp -s "$a" "$b"; then
        echo "- [FAIL] ${label}: files differ (${a} vs ${b})"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

assert_file_missing() {
    local label="$1" f="$2"
    if [[ -e "$f" ]]; then
        echo "- [FAIL] ${label}: ${f} should not exist yet"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

echo ">> Test incremental finalize (§6.5/§11, per-method augmented)"

tmpdir=$(mktemp -d)
export CFG_OUTDIR="${tmpdir}/out"
export CFG_CHROMOSOMES="22"
export CFG_METHODS="sbayesr ldpred2"

sumstat_name="trait_smoke"
sumstat_dir="${CFG_OUTDIR}/sumstats/${sumstat_name}"
prep_dir="${CFG_OUTDIR}/prep"
work="${sumstat_dir}/work"
mkdir -p "$prep_dir" "$work"

# Prep-style map: col4=geno_snpid, col7=ldref_snpid (finalize vm_base join key)
cat > "${sumstat_dir}/variant_map.tsv" <<'EOF'
chr	pos_b37	pos_b38	geno_snpid	geno_a1	geno_a2	ldref_snpid	ldref_a1	ldref_a2	ldref_a2freq
22	17309296	16828406	22:16828406_A_G	A	G	rs175146	G	A	0.11
22	18181984	17699218	22:17699218_T_G	T	G	rs17207051	G	T	0.02
EOF

# LDpred2 variant map (12-col prep schema) shares the same rsids here.
cat > "${prep_dir}/variant_map_ldpred2.tsv" <<'EOF'
chr	pos_b37	pos_b38	geno_snpid	geno_a1	geno_a2	ldref_snpid	ldref_a1	ldref_a2	ldref_a2freq	ld	block_id
22	17309296	16828406	22:16828406_A_G	A	G	rs175146	G	A	0.11	0.1	1
22	18181984	17699218	22:17699218_T_G	T	G	rs17207051	G	T	0.02	0.1	1
EOF

mkdir -p "${work}/filtered_sbayesr" "${work}/filtered_ldpred2"
for d in filtered_sbayesr filtered_ldpred2; do
cat > "${work}/${d}/chr22_matched.tsv" <<'EOF'
LDREF_SNPID	B	SE	Z	P
rs175146	-0.01	0.02	-0.5	0.6
rs17207051	0.02	0.03	0.67	0.5
EOF
done

mkdir -p "${work}/posteriors_mapped"
cat > "${work}/posteriors_mapped/chr22.snpRes" <<'EOF'
ID	A1	A2	Freq	Effect	SE	PIP
22:16828406_A_G	A	G	0.11	-0.001	0.005	0.61
22:17699218_T_G	T	G	0.02	-0.002	0.008	0.65
EOF

mkdir -p "${work}/scores"
cat > "${work}/scores/chr22.sscore" <<'EOF'
#FID	IID	SCORE1_SUM	ALLELE_CT	N_VARIANTS
0	sample1	1.23	100	2
EOF

# Day 1: sBayesR only
combine_scores_for_method "$sumstat_name" "sbayesr"
cp "${sumstat_dir}/scores_sbayesr.gz" "${tmpdir}/day1_scores_sbayesr.gz"

write_augmented_sumstat_for_method "$sumstat_dir" "$prep_dir" "sbayesr"
day1_header=$(zcat "${sumstat_dir}/augmented_sbayesr.gz" | head -1)
assert_contains "day1 sbayesr postEffect" "$day1_header" "postEffect"
assert_contains "day1 sbayesr benchEffect" "$day1_header" "benchEffect"
assert_not_contains "day1 sbayesr no postp" "$day1_header" "postp_ldpred2"
assert_file_missing "day1 no augmented_ldpred2.gz" "${sumstat_dir}/augmented_ldpred2.gz"
cp "${sumstat_dir}/augmented_sbayesr.gz" "${tmpdir}/day1_augmented_sbayesr.gz"

# Day 2: add LDpred2 posteriors + scores; re-combine only ldpred2
mkdir -p "${work}/posteriors_mapped_ldpred2" "${work}/scores_ldpred2"
cat > "${work}/posteriors_mapped_ldpred2/chr22.snpRes" <<'EOF'
ID	A1	A2	Freq	Effect	SE	PIP
22:16828406_A_G	A	G	0.11	-0.0015	0.004	0.70
22:17699218_T_G	T	G	0.02	-0.0025	0.007	0.55
EOF
cat > "${work}/scores_ldpred2/chr22.sscore" <<'EOF'
#FID	IID	SCORE1_SUM	ALLELE_CT	N_VARIANTS
0	sample1	1.45	100	2
EOF

combine_scores_for_method "$sumstat_name" "ldpred2"
assert_eq_file "scores_sbayesr.gz unchanged" \
    "${tmpdir}/day1_scores_sbayesr.gz" "${sumstat_dir}/scores_sbayesr.gz"

write_augmented_sumstat_for_method "$sumstat_dir" "$prep_dir" "ldpred2"
day2_header=$(zcat "${sumstat_dir}/augmented_ldpred2.gz" | head -1)
assert_contains "day2 ldpred2 postEffect" "$day2_header" "postEffect"
assert_contains "day2 ldpred2 postp_ldpred2" "$day2_header" "postp_ldpred2"
assert_contains "day2 ldpred2 benchEffect" "$day2_header" "benchEffect"
assert_eq_file "augmented_sbayesr.gz unchanged" \
    "${tmpdir}/day1_augmented_sbayesr.gz" "${sumstat_dir}/augmented_sbayesr.gz"

if [[ ! -f "${sumstat_dir}/scores_ldpred2.gz" ]]; then
    echo "- [FAIL] scores_ldpred2.gz missing after day2 combine"
    failures=$((failures + 1))
else
    echo "- [OK] scores_ldpred2.gz created"
fi

if [[ $failures -gt 0 ]]; then
    echo "Tests failed: ${failures}"
    exit 1
fi
echo "All incremental finalize tests passed."
