#!/usr/bin/env bash

set -euo pipefail

test_script="find_col_indices"
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

  "${test_script}.sh" ./sumstat "RSID,EffectAllele,B" > ./observed-results.tsv

  _check_results ./observed-results.tsv ./expected-result_1.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# base case sbayesr
_setup "base_case"

##### REUSE THE SCRIPT YOU JUT WROTE it should work!!
cat <<EOF > ./sumstat
CHR	POS	CHR	POS	0	RSID	EffectAllele	OtherAllele	B	SE	Z	P	OR	Neff	CaseN	ControlN	EAF	INFO	N	
2	72560667	2	72333538	4630538	rs75240917	C	T	0.076896	0.026	2.957048	0.003106	1.07993	121702.873174	52301	72744	0.97	0.868	111487
2	9929213	2	9789084	3044895	rs10929615	T	G	0.008196	0.0092	0.898160	0.3691	1.00823	150092.970769	64952	88856	0.71	0.901	111487
2	174849052	2	173984324	876359	rs186552910	C	T	0.050703	0.0435	1.165047	0.244	1.05201	97308.179613	41354	59084	0.99	0.923	111487
2	102000790	2	101384328	7273268	rs62156072	G	A	0.005803	0.0158	0.364747	0.7153	1.00582	150092.970769	64952	88856	0.95	0.96	111487
2	195451582	2	194586858	4104750	rs17587775	T	C	-0.006904	0.021	-0.331191	0.7405	0.99312	150092.970769	64952	88856	0.97	0.979	111487
EOF

cat <<EOF > ./expected-result_1.tsv
6,7,9
EOF

_run_script

