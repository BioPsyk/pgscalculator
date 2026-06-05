#!/usr/bin/env bash
# Unit test for the LDpred2 variant-map build (prep-inclusion-list-ldpred2).
#
# Covers the parallelized layout introduced in v2.2.x:
#   - single-chromosome (SLURM array) mode: writes a per-chr partial + marker,
#     does NOT produce the combined map;
#   - the combine step concatenates partials into prep/variant_map_ldpred2.tsv;
#   - the full-run path (per-chr partials in parallel + auto-combine) produces
#     a byte-identical combined map;
#   - join semantics: direct / allele-swap / strand-complement matches, NA on
#     genotype miss, af_UKBB -> ldref_a2freq reconciliation (1 - af), and the
#     ld / block_id columns carried through;
#   - opt-in: with ldpred2.ld_dir unset the step is a no-op.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
export SCRIPT_DIR="${PROJECT_DIR}/bin"
export LIB_DIR="${PROJECT_DIR}/bin/lib"
export STEPS_DIR="${LIB_DIR}/steps"
VERBOSE=0
FORCE=0

source "${LIB_DIR}/common.sh"
source "${STEPS_DIR}/prep_inclusion_list_ldpred2.sh"

tmpdir=""
failures=0

cleanup() { [[ -n "$tmpdir" && -d "$tmpdir" ]] && rm -rf "$tmpdir"; }
trap cleanup EXIT

assert_eq_file() {
    local label="$1" obs="$2" exp="$3"
    if ! diff -u "$exp" "$obs" > "${tmpdir}/diff.out" 2>&1; then
        echo "- [FAIL] ${label}: output differs from expected"
        cat "${tmpdir}/diff.out"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

assert_true() {
    local label="$1"; shift
    if "$@"; then echo "- [OK] ${label}"; else echo "- [FAIL] ${label}"; failures=$((failures + 1)); fi
}

assert_false() {
    local label="$1"; shift
    if "$@"; then echo "- [FAIL] ${label}"; failures=$((failures + 1)); else echo "- [OK] ${label}"; fi
}

echo ">> Test variant_map_for_ldpred2 (prep-inclusion-list-ldpred2)"

tmpdir=$(mktemp -d)

# ---- shared fixtures -------------------------------------------------------
# Build genotypes + ldref_ldpred2/map.tsv under a prep dir for a given outdir.
make_fixtures() {
    local outdir="$1"
    local prep="${outdir}/prep"
    mkdir -p "${prep}/genotypes" "${prep}/ldref_ldpred2"

    # snplist_sorted only needs to exist (deps check); content irrelevant here.
    echo "rs1" > "${prep}/genotypes/snplist_sorted"

    # chr22_pvar_fmt: chr:pos <TAB> a1 <TAB> a2 <TAB> snpid  (GRCh37 positions)
    # v1 direct (G/A), v2 swap (stored C/T vs ldref T/C), v3 complement
    # (stored C/T vs ldref G/A), v4 absent -> NA in map.
    {
        printf '22:100\tG\tA\tsnp_v1\n'
        printf '22:300\tC\tT\tsnp_v2\n'
        printf '22:500\tC\tT\tsnp_v3\n'
    } > "${prep}/genotypes/chr22_pvar_fmt"

    # map.tsv: chr pos_b37 pos_b38 a0 a1 rsid af_UKBB ld block_id
    # a1 = effect/ALT, a0 = other/REF; ldref_a2freq = 1 - af_UKBB.
    {
        printf 'chr\tpos_b37\tpos_b38\ta0\ta1\trsid\taf_UKBB\tld\tblock_id\n'
        printf '22\t100\t200\tA\tG\trs1\t0.30\t1.5\t5\n'
        printf '22\t300\t400\tC\tT\trs2\t0.10\t2.0\t5\n'
        printf '22\t500\t600\tA\tG\trs3\t0.25\t1.1\t6\n'
        printf '22\t700\t800\tA\tT\trs4\t0.05\t0.9\t6\n'
    } > "${prep}/ldref_ldpred2/map.tsv"
}

# Expected per-chr partial (no header). 12 cols:
# chr pos_b37 pos_b38 geno_snpid geno_a1 geno_a2 ldref_rsid ldref_a1 ldref_a2 ldref_a2freq ld block_id
expected_partial="${tmpdir}/expected_chr22.map"
{
    printf '22\t100\t200\tsnp_v1\tG\tA\trs1\tG\tA\t0.7\t1.5\t5\n'
    printf '22\t300\t400\tsnp_v2\tC\tT\trs2\tT\tC\t0.9\t2.0\t5\n'
    printf '22\t500\t600\tsnp_v3\tG\tA\trs3\tG\tA\t0.75\t1.1\t6\n'
    printf '22\t700\t800\tNA\tNA\tNA\trs4\tT\tA\t0.95\t0.9\t6\n'
} > "$expected_partial"

# Expected combined map (header + body).
expected_combined="${tmpdir}/expected_combined.tsv"
{
    printf 'chr\tpos_b37\tpos_b38\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq\tld\tblock_id\n'
    cat "$expected_partial"
} > "$expected_combined"

export CFG_GENOTYPE_BUILD="GRCh37"
export CFG_LDPRED2_LD_DIR="dummy"   # opt-in flag (presence only)

# ===========================================================================
# Scenario A: per-chromosome array task (single_chr_task=1) -> partial + combine
# ===========================================================================
echo "-- scenario A: array task (single_chr_task=1) + combine"
outA="${tmpdir}/outA"
make_fixtures "$outA"
export CFG_OUTDIR="$outA"
export CFG_CHROMOSOMES="22"
export CFG_SINGLE_CHR_TASK="1"       # driver array task -> partial-only

run_prep_inclusion_list_ldpred2 > "${tmpdir}/A.log" 2>&1 || { cat "${tmpdir}/A.log"; echo "run failed"; exit 1; }

assert_eq_file "A: per-chr partial matches expected" \
    "${outA}/prep/inclusion_list_ldpred2/chr22.map" "$expected_partial"
assert_true  "A: per-chr completion marker written" \
    test -f "${outA}/prep/inclusion_list_ldpred2/.completed_chr22"
assert_false "A: combined map NOT created in single-chr mode" \
    test -f "${outA}/prep/variant_map_ldpred2.tsv"

run_prep_inclusion_list_ldpred2_combine > "${tmpdir}/A_combine.log" 2>&1 || { cat "${tmpdir}/A_combine.log"; echo "combine failed"; exit 1; }

assert_eq_file "A: combined map matches expected" \
    "${outA}/prep/variant_map_ldpred2.tsv" "$expected_combined"
assert_true  "A: combine wrote step .completed" \
    test -f "${outA}/prep/inclusion_list_ldpred2/.completed"
assert_true  "A: combine wrote details tsv" \
    test -f "${outA}/prep/details/prep_inclusion_list_ldpred2_steps.tsv"

# details should report 3 matched of 4 (75.00%)
if grep -q "match_rate=75.00%" "${outA}/prep/details/prep_inclusion_list_ldpred2_steps.tsv"; then
    echo "- [OK] A: details report 75.00% match rate"
else
    echo "- [FAIL] A: details match rate"; cat "${outA}/prep/details/prep_inclusion_list_ldpred2_steps.tsv"; failures=$((failures + 1))
fi

# ===========================================================================
# Scenario B: single-chromosome CONFIG without the array flag must run the full
# prep (incl. combine). Regression guard: a `chromosomes: 22` config in a normal
# local run previously fell into partial-only mode and never combined.
# ===========================================================================
echo "-- scenario B: single-chr config (no flag) runs full prep + combine"
outB="${tmpdir}/outB"
make_fixtures "$outB"
export CFG_OUTDIR="$outB"
export CFG_CHROMOSOMES="22"           # single value...
unset CFG_SINGLE_CHR_TASK             # ...but NOT an array task -> full run

run_prep_inclusion_list_ldpred2 > "${tmpdir}/B.log" 2>&1 || { cat "${tmpdir}/B.log"; echo "run failed"; exit 1; }

assert_eq_file "B: full-run combined map matches expected" \
    "${outB}/prep/variant_map_ldpred2.tsv" "$expected_combined"
assert_true  "B: full-run wrote step .completed" \
    test -f "${outB}/prep/inclusion_list_ldpred2/.completed"

# A and B combined maps must be byte-identical.
assert_eq_file "B: full-run == array+combine" \
    "${outB}/prep/variant_map_ldpred2.tsv" "${outA}/prep/variant_map_ldpred2.tsv"
export CFG_SINGLE_CHR_TASK="1"

# ===========================================================================
# Scenario C: opt-in skip when ldpred2.ld_dir is unset
# ===========================================================================
echo "-- scenario C: opt-in skip"
outC="${tmpdir}/outC"
make_fixtures "$outC"
export CFG_OUTDIR="$outC"
export CFG_CHROMOSOMES="22"
unset CFG_LDPRED2_LD_DIR

run_prep_inclusion_list_ldpred2 > "${tmpdir}/C.log" 2>&1 || { cat "${tmpdir}/C.log"; echo "run failed"; exit 1; }
assert_false "C: no partial written when ld_dir unset" \
    test -f "${outC}/prep/inclusion_list_ldpred2/chr22.map"
assert_false "C: no combined map when ld_dir unset" \
    test -f "${outC}/prep/variant_map_ldpred2.tsv"
export CFG_LDPRED2_LD_DIR="dummy"

if [[ $failures -gt 0 ]]; then
    echo "Tests failed: ${failures}"
    exit 1
fi
echo "All variant_map_for_ldpred2 tests passed."
