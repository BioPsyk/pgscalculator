#!/usr/bin/env bash

set -euo pipefail

test_script="make_geno_pvar_snpid_unique"
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
  "${test_script}.sh" ./input.pvar ./observed-results.pvar
  _check_results ./observed-results.pvar ./expected-result_1.pvar

  echo "- [OK] ${curr_case}"
  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Base case - Make SNP IDs unique

_setup "base_case"

cat <<-EOF > ./input.pvar
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
##contig=<ID=22,assembly=b37,length=51304566>
##ALT=<ID=SNP,Description="SNP">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16050075	rs123	A	G	.	PASS	.
22	16050115	rs123	G	A	.	PASS	.
22	16050213	rs789	C	T	.	PASS	.
22	16050319	.	C	T	.	PASS	.
22	16050527	rs321	C	A	.	PASS	.
EOF

cat <<-EOF > ./expected-result_1.pvar
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
##contig=<ID=22,assembly=b37,length=51304566>
##ALT=<ID=SNP,Description="SNP">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16050075	rs123_1	A	G	.	PASS	.
22	16050115	rs123_2	G	A	.	PASS	.
22	16050213	rs789	C	T	.	PASS	.
22	16050319	.	C	T	.	PASS	.
22	16050527	rs321	C	A	.	PASS	.
EOF

_run_script 