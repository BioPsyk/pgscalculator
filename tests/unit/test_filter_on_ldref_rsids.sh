#!/usr/bin/env bash

set -euo pipefail

test_script="filter_on_ldref_rsids"
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

  "${test_script}.sh" ./input_file.tsv.gz ./rsids.txt > ./observed-results.tsv
  _check_results ./observed-results.tsv ./expected-result_1.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# keep rsids

_setup "keep two rsid rows"

cat <<EOF > ./rsids.txt
rs6831643
rs12726220
EOF

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
0	RSID	CHR	POS	EffectAllele	OtherAllele	EAF	B	SE	P
100	rs6831643	4	99833465	T	C	0.669	-0.0321	0.0137	0.0193
1003	rs12726220	1	150984623	A	G	0.948	-0.0315	0.0277	0.2547
EOF

_run_script

