#!/usr/bin/env bash

set -euo pipefail

test_script="prepare_benchmark_scoring"
initial_dir=$(pwd)"/${test_script}"
curr_case=""

mkdir "${initial_dir}"
cd "${initial_dir}"

#=================================================================================
# Helpers
#=================================================================================

function _setup {
  mkdir "${1}"
  cd "${1}"
  curr_case="${1}"
}

function _run_script {
  bash ../../bin/prepare_benchmark_scoring.sh "${@}" > "./result_1.txt"
  
  if ! diff "result_1.txt" "expected-result_1.tsv"; then
    echo "Failed: ${curr_case}"
    exit 1
  fi
}

#---------------------------------------------------------------------------------
# Case 1: Closely positioned variants (should prune)
_setup "close_variants"

# Create input sumstat file with closely positioned variants
cat <<EOF > ./sumstat
CHR POS RSID EffectAllele OtherAllele B SE P
20 1000000 rs1 A G 0.1 0.02 0.001
20 1000100 rs2 T C -0.05 0.01 0.002
20 1000200 rs3 C G 0.15 0.03 0.003
20 1000300 rs4 A T -0.08 0.02 0.004
20 1000400 rs5 G C 0.12 0.03 0.005
20 1000500 rs6 T A -0.06 0.02 0.006
20 1000600 rs7 C T 0.09 0.02 0.007
20 1000700 rs8 G A -0.11 0.03 0.008
20 1000800 rs9 A C 0.07 0.02 0.009
20 1000900 rs10 T G -0.14 0.03 0.010
EOF

# Create PLINK files
cat <<EOF > ./geno.pvar
#CHROM  POS     ID      REF     ALT
20      1000000 rs1     G       A
20      1000100 rs2     C       T
20      1000200 rs3     G       C
20      1000300 rs4     T       A
20      1000400 rs5     C       G
20      1000500 rs6     A       T
20      1000600 rs7     T       C
20      1000700 rs8     A       G
20      1000800 rs9     C       A
20      1000900 rs10    G       T
EOF

cat <<EOF > ./geno.psam
#IID    SEX
sample1 1
sample2 2
EOF

# Expected output after pruning (subset of variants)
cat <<EOF > ./expected-result_1.tsv
CHR POS RSID EffectAllele OtherAllele B SE P
20 1000000 rs1 A G 0.1 0.02 0.001
20 1000400 rs5 G C 0.12 0.03 0.005
20 1000800 rs9 A C 0.07 0.02 0.009
EOF

_run_script sumstat geno 20

#---------------------------------------------------------------------------------
# Case 2: Distant variants (should reproduce error)
_setup "distant_variants"

# Create input sumstat file with distant variants
cat <<EOF > ./sumstat
CHR POS RSID EffectAllele OtherAllele B SE P
20 1000000 rs1 A G 0.1 0.02 0.001
20 20000000 rs2 T C -0.05 0.01 0.002
20 40000000 rs3 C G 0.15 0.03 0.003
20 60000000 rs4 A T -0.08 0.02 0.004
20 80000000 rs5 G C 0.12 0.03 0.005
EOF

# Create PLINK files
cat <<EOF > ./geno.pvar
#CHROM  POS     ID      REF     ALT
20      1000000 rs1     G       A
20      20000000 rs2    C       T
20      40000000 rs3    G       C
20      60000000 rs4    T       A
20      80000000 rs5    C       G
EOF

cat <<EOF > ./geno.psam
#IID    SEX
sample1 1
sample2 2
EOF

# This case should reproduce the error we're seeing in production
# We expect PLINK to fail the pruning due to no variants in LD

_run_script sumstat geno 20 || echo "Expected error occurred in distant_variants test"
