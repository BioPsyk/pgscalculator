#!/usr/bin/env bash

set -euo pipefail

test_script="concatenate_plink_maf"
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
  "${test_script}.sh" "raw_maf_chrall" ./chr*.afreq
  _check_results "raw_maf_chrall" ./expected-result_1.txt

  echo "- [OK] ${curr_case}"
  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Standard PLINK2 .afreq format

_setup "standard_plink2_format"

cat <<EOF > ./chr1.afreq
#CHROM	ID	REF	ALT	PROVISIONAL_REF?	ALT_FREQS	OBS_CT
1	rs123456	A	T	Y	0.25	1000
1	rs789012	C	G	Y	0.35	1000
EOF

cat <<EOF > ./chr2.afreq
#CHROM	ID	REF	ALT	PROVISIONAL_REF?	ALT_FREQS	OBS_CT
2	rs345678	G	A	Y	0.15	1000
2	rs901234	T	C	Y	0.45	1000
EOF

cat <<EOF > ./expected-result_1.txt
CHR SNP A1 A2 MAF NCHROBS
1 rs123456 A T 0.25 1000
1 rs789012 C G 0.35 1000
2 rs345678 G A 0.15 1000
2 rs901234 T C 0.45 1000
EOF

_run_script

#---------------------------------------------------------------------------------
# Different column order

_setup "different_column_order"

cat <<EOF > ./chr3.afreq
#ID	CHROM	ALT_FREQS	REF	ALT	OBS_CT	EXTRA_COL
rs111111	3	0.55	A	G	1000	X
rs222222	3	0.65	T	A	1000	Y
EOF

cat <<EOF > ./expected-result_1.txt
CHR SNP A1 A2 MAF NCHROBS
3 rs111111 A G 0.55 1000
3 rs222222 T A 0.65 1000
EOF

_run_script

#---------------------------------------------------------------------------------
# Alternative CHR column name

_setup "alternative_chr_column"

cat <<EOF > ./chr4.afreq
#CHR	ID	REF	ALT	ALT_FREQS	OBS_CT
4	rs333333	C	T	0.75	1000
4	rs444444	G	C	0.85	1000
EOF

cat <<EOF > ./expected-result_1.txt
CHR SNP A1 A2 MAF NCHROBS
4 rs333333 C T 0.75 1000
4 rs444444 G C 0.85 1000
EOF

_run_script
