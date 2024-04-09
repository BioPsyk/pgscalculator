#!/usr/bin/env bash

set -euo pipefail

test_script="sumstat_maf_filter"
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
  "${test_script}.sh" "./ssfile" "./maffile" "0.05" > "observed" 

  _check_results observed ./expected-result_1.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"
#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Reformat_prscs

_setup "base case maf"

cat <<EOF > ./ssfile
posteriorID	EA	Effect	genoID
rs175146	A	0.006499	rs175146
rs17207051	T	-0.016404	rs17207051
EOF

cat <<EOF > ./maffile
CHR SNP A1 A2 MAF NCHROBS
22 rs175146 G A 0.1711 5008
22 rs17207051 T C 0.09026 5008
22 rs71328205 C T 0.03431 5008
EOF

cat <<EOF > ./expected-result_1.tsv
posteriorID	EA	Effect	genoID
rs175146	A	0.006499	rs175146
rs17207051	T	-0.016404	rs17207051
EOF

_run_script


