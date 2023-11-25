#!/usr/bin/env bash

set -euo pipefail

test_script="filter_bad_values"
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

  "${test_script}.sh" ./input_file.tsv.gz > ./observed-results.tsv
  _check_results ./observed-results.tsv ./expected-result_1.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Add n

_setup "add beta and se, if missing"

cat <<EOF | gzip -c > ./input_file.tsv.gz
0	RSID	CHR	POS	EffectAllele	OtherAllele	EAF	P	N	Z	B	SE
1	rs6439928	3	141663261	T	C	0.658	0.2648	13688	0.15	0.00191108	0.0127406
10	rs6463169	7	42980893	T	C	0.825	0.2012	13688	0.22	0.00349938	0.0159063
100	rs6831643	4	99833465	T	C	0.0	0.0193	13688	0.11	0.0014128	0.0128436
1000	rs10197378	2	29092758	A	G	0.183	0.2226	13688	-0.24	0	0.0156307
1001	rs10021082	4	100801356	T	C	0.958	0.2265	13688	0.12	0.00361567	0
1002	rs12709653	18	27735538	A	G	0.775	0.3176	13688	0.23	0.00332889	0.0144734
1003	rs12726220	1	150984623	A	G	1.00	0.2547	13688	-0.12	-0.00326656	0.0272213
1004	rs12739293	1	118812591	T	C	0.133	0.6873	13688	-0.21	-0.00373765	0.0177983
1005	rs12754538	1	8408079	T	C	0.308	-6e-04	13688	-0.05	-0.000654571	0.0130914
EOF

cat <<EOF > ./expected-result_1.tsv
0	RSID	CHR	POS	EffectAllele	OtherAllele	EAF	P	N	Z	B	SE
1	rs6439928	3	141663261	T	C	0.658	0.2648	13688	0.15	0.00191108	0.0127406
10	rs6463169	7	42980893	T	C	0.825	0.2012	13688	0.22	0.00349938	0.0159063
1002	rs12709653	18	27735538	A	G	0.775	0.3176	13688	0.23	0.00332889	0.0144734
1004	rs12739293	1	118812591	T	C	0.133	0.6873	13688	-0.21	-0.00373765	0.0177983
1005	rs12754538	1	8408079	T	C	0.308	-6e-04	13688	-0.05	-0.000654571	0.0130914
EOF

_run_script

