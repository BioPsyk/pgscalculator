#!/usr/bin/env bash

set -euo pipefail

test_script="make_augmented_output"
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

  "${test_script}.sh" ./sumstat ./map ./maffile ./posterior ./benchmark > "observed.tsv" 

  _check_results observed.tsv ./expected-result_1.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Reformat_prscs

_setup "base case"

cat <<EOF > ./sumstat
CHR	POS	0	RSID	EffectAllele	OtherAllele	EAF	B	SE	P	Z
22	141663261	1	rs5992124	T	C	0.658	-0.0157	0.0141	0.2648	-0.0141
22	100801356	2	rs17207051	T	C	0.958	0.0319	0.0264	0.2265	0.0264
22	27735538	3	rs78872431	A	G	0.775	-0.0142	0.0142	0.3176	-0.0142
22	150984623	4	rs8138892	A	G	0.948	-0.0315	0.0277	0.2547	-0.0277
22	8408079	5	rs71328205	T	C	0.308	-6e-04	0.015	0.9663	-0.9663
EOF

cat <<EOF > ./map
b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
22:17434931	22:16954041	rs8138892	C	T	rs8138892	T	C	NA	NA	NA
22:17685648	22:17204758	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
22:17697944	22:17217054	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
22:18181984	22:17699218	rs5992124	T	G	rs175146	G	T	rs17207051	G	T
22:18334573	22:17851807	rs17207051	T	C	rs17207051	C	T	rs5992124	C	T
EOF

cat <<EOF > ./maffile
CHR SNP A1 A2 MAF NCHROBS
22 rs175146 G A 0.1711 5008
22 rs17207051 T C 0.09026 5008
22 rs71328205 C T 0.03431 5008
22 rs78872431 C T 0.0531 5008
22 rs8138892 C T 0.04531 5008
EOF

cat <<EOF > ./posterior
posteriorID	EA	Effect	genoID
rs175146	A	0.006499	rs175146
rs17207051	T	-0.016404	rs17207051
EOF

cat <<EOF > ./benchmark
posteriorID	EA	Effect	genoID
rs175146	A	0.005399	rs175146
rs17207051	T	-0.019404	rs17207051
EOF


# SNP     A1      A2      BETA    P
cat <<EOF > ./expected-result_1.tsv
RSID	CHR	POS	EffectAllele	OtherAllele	B	SE	Z	P	genoID	MAF	postEffect	benchEffect
rs5992124	22	141663261	T	C	-0.0157	0.0141	-0.0141	0.2648	rs175146	0.1711	0.006499	0.005399
rs17207051	22	100801356	T	C	0.0319	0.0264	0.0264	0.2265	rs17207051	0.09026	-0.016404	-0.019404
rs78872431	22	27735538	A	G	-0.0142	0.0142	-0.0142	0.3176	rs78872431	0.0531	NA	NA
rs8138892	22	150984623	A	G	-0.0315	0.0277	-0.0277	0.2547	rs8138892	0.04531	NA	NA
rs71328205	22	8408079	T	C	-6e-04	0.015	-0.9663	0.9663	rs71328205	0.03431	NA	NA
EOF

_run_script


