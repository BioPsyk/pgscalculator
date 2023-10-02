#!/usr/bin/env bash

set -euo pipefail

test_script="format_sumstats"
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

mapfile="${ch_assets_sumstats_header_map}"
function _run_script {
  method=$1

  "${test_script}.sh" input_file.tsv ${mapfile} ${method} "formatted" 

  _check_results formatted_1.tsv ./expected-result_1.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Reformat_prscs

_setup "Reformat prscs base case"

cat <<EOF > ./input_file.tsv
0	RSID	CHR	POS	EffectAllele	OtherAllele	EAF	B	SE	P
1	rs6439928	3	141663261	T	C	0.658	-0.0157	0.0141	0.2648
10	rs6463169	7	42980893	T	C	0.825	-0.0219	0.0171	0.2012
100	rs6831643	4	99833465	T	C	0.669	-0.0321	0.0137	0.0193
1000	rs10197378	2	29092758	A	G	0.183	-0.0189	0.0155	0.2226
1001	rs10021082	4	100801356	T	C	0.958	0.0319	0.0264	0.2265
1002	rs12709653	18	27735538	A	G	0.775	-0.0142	0.0142	0.3176
1003	rs12726220	1	150984623	A	G	0.948	-0.0315	0.0277	0.2547
1004	rs12739293	1	118812591	T	C	0.133	0.007	0.0175	0.6873
1005	rs12754538	1	8408079	T	C	0.308	-6e-04	0.015	0.9663
EOF

cat <<EOF > ./expected-result_1.tsv
SNP	A1	A2	BETA	P
rs12726220	A	G	-0.0315	0.2547
rs12739293	T	C	0.007	0.6873
rs12754538	T	C	-6e-04	0.9663
EOF

_run_script "prscs"

#---------------------------------------------------------------------------------
# Next case

#_setup "valid_rows_missing_afreq"
#
#cat <<EOF > ./acor.tsv
#0	A1	A2	CHRPOS	RSID	EffectAllele	OtherAllele	EMOD
#1	A	G	12:126406434	rs1000000	G	A	-1
#EOF
#
#cat <<EOF > ./stat.tsv
#0	B	SE	Z	P
#1	-0.0143	0.0156	-0.916667	0.3604
#EOF
#
#_run_script
