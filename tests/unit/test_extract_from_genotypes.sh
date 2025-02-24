#!/usr/bin/env bash

set -euo pipefail

test_script="extract_from_genotypes"
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

function _check_results {
  obs=$1
  exp=$2
  if ! diff ${obs} ${exp} &> ./difference; then
    echo "- [FAIL] ${curr_case}"
    cat ./difference 
    exit 1
  fi
}

function _run_script {
  "${test_script}.sh" ./input.pgen ./input.pvar ./input.psam ./extract.txt ./output

  # Check PVAR file
  _check_results ./output.pvar ./expected.pvar
  # Check PSAM file
  _check_results ./output.psam ./expected.psam

  echo "- [OK] ${curr_case}"
  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Base case - Extract subset of variants

_setup "base_case"

# Create a dummy binary PGEN file (we can't easily create real PGEN content in a test)
dd if=/dev/zero of=./input.pgen bs=1 count=100

cat <<EOF > ./input.pvar
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
##contig=<ID=22,assembly=b37,length=51304566>
##ALT=<ID=SNP,Description="SNP">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
22      16050075        rs123   A       G       .       PASS    .
22      16050115        rs456   G       A       .       PASS    .
22      16050213        rs789   C       T       .       PASS    .
22      16050319        rs101   C       T       .       PASS    .
22      16050527        rs321   C       A       .       PASS    .
EOF

cat <<EOF > ./input.psam
#FID    IID     SEX
fam1    ind1    1
fam1    ind2    2
fam2    ind3    2
fam2    ind4    1
EOF

# Create extract file with subset of variants to keep
cat <<EOF > ./extract.txt
rs123
rs789
rs321
EOF

# Expected PVAR output after extraction
cat <<EOF > ./expected.pvar
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
22      16050075        rs123   A       G       .       PASS    .
22      16050213        rs789   C       T       .       PASS    .
22      16050527        rs321   C       A       .       PASS    .
EOF

# Expected PSAM output (should be unchanged)
cat <<EOF > ./expected.psam
#FID    IID     SEX
fam1    ind1    1
fam1    ind2    2
fam2    ind3    2
fam2    ind4    1
EOF

_run_script

#---------------------------------------------------------------------------------
# Case with missing variants

_setup "missing_variants"

# Create a dummy binary PGEN file
dd if=/dev/zero of=./input.pgen bs=1 count=100

cat <<EOF > ./input.pvar
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
##contig=<ID=22,assembly=b37,length=51304566>
##ALT=<ID=SNP,Description="SNP">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
22      16050075        rs123   A       G       .       PASS    .
22      16050115        rs456   G       A       .       PASS    .
22      16050213        rs789   C       T       .       PASS    .
EOF

cat <<EOF > ./input.psam
#FID    IID     SEX
fam1    ind1    1
fam1    ind2    2
EOF

# Create extract file with some variants that don't exist
cat <<EOF > ./extract.txt
rs123
rs999
rs789
rs888
EOF

# Expected PVAR output after extraction
cat <<EOF > ./expected.pvar
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
22      16050075        rs123   A       G       .       PASS    .
22      16050213        rs789   C       T       .       PASS    .
EOF

# Expected PSAM output (should be unchanged)
cat <<EOF > ./expected.psam
#FID    IID     SEX
fam1    ind1    1
fam1    ind2    2
EOF

_run_script 