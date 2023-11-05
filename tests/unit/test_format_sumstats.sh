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

  "${test_script}.sh" input_file.tsv.gz ${mapfile} ${method} "formatted" 
  _check_results formatted_1.tsv ./expected-result_1.tsv
  _check_results formatted_2.tsv ./expected-result_2.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Reformat_prscs

_setup "prscs base case"

cat <<EOF | gzip -c > ./input_file.tsv.gz
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

cat <<EOF > ./expected-result_2.tsv
SNP	A1	A2	BETA	P
rs10197378	A	G	-0.0189	0.2226
EOF

_run_script "prscs"

#---------------------------------------------------------------------------------
# Next case

_setup "Reformat sbayesr base case"

cat <<EOF | gzip -c > ./input_file.tsv.gz
0	RSID	CHR	POS	EffectAllele	OtherAllele	EAF	B	SE	P	N
1	rs6439928	3	141663261	T	C	0.658	-0.0157	0.0141	0.2648	30000
10	rs6463169	7	42980893	T	C	0.825	-0.0219	0.0171	0.2012	30000
100	rs6831643	4	99833465	T	C	0.669	-0.0321	0.0137	0.0193	30000
1000	rs10197378	2	29092758	A	G	0.183	-0.0189	0.0155	0.2226	30000
1001	rs10021082	4	100801356	T	C	0.958	0.0319	0.0264	0.2265	30000
1002	rs12709653	18	27735538	A	G	0.775	-0.0142	0.0142	0.3176	30000
1003	rs12726220	1	150984623	A	G	0.948	-0.0315	0.0277	0.2547	30000
1004	rs12739293	1	118812591	T	C	0.133	0.007	0.0175	0.6873	30000
1005	rs12754538	1	8408079	T	C	0.308	-6e-04	0.015	0.9663	30000
EOF

cat <<EOF > ./expected-result_1.tsv
SNP A1 A2 freq b se p N
rs12726220 A G 0.948 -0.0315 0.0277 0.2547 30000
rs12739293 T C 0.133 0.007 0.0175 0.6873 30000
rs12754538 T C 0.308 -6e-04 0.015 0.9663 30000
EOF

cat <<EOF > ./expected-result_2.tsv
SNP A1 A2 freq b se p N
rs10197378 A G 0.183 -0.0189 0.0155 0.2226 30000
EOF

_run_script "sbayesr"
