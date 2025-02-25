#!/usr/bin/env bash

set -euo pipefail

test_script="make_snplist_from_pvar"
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
  "${test_script}.sh" ./input*.pvar > ./observed-results.txt
  _check_results ./observed-results.txt ./expected-result_1.txt

  echo "- [OK] ${curr_case}"
  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Base case - Extract SNP IDs from multiple PVAR files

_setup "base_case"

cat <<EOF > ./input1.pvar
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
##contig=<ID=22,assembly=b37,length=51304566>
##ALT=<ID=SNP,Description="SNP">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="EAS Allele Frequency">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="EUR Allele Frequency">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="AFR Allele Frequency">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="AMR Allele Frequency">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="SAS Allele Frequency">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
22      16050075        rs123   A       G       .       PASS    .
22      16050115        rs456   G       A       .       PASS    .
22      16050213        rs789   C       T       .       PASS    .
EOF

cat <<EOF > ./input2.pvar
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
22      16050319        rs101   C       T       .       PASS    .
22      16050527        rs321   C       A       .       PASS    .
EOF

cat <<EOF > ./expected-result_1.txt
rs101
rs123
rs321
rs456
rs789
EOF

_run_script 