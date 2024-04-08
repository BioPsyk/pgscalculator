#!/usr/bin/env bash

set -euo pipefail

test_script="format_posteriors"
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

  "${test_script}.sh" ./posteriors "$1" ./mapfile_noNA "$2" "$3" > ./observed-results.tsv

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

snp_posteriors_cols="2,5,8"
map_from_to="3,6"

cat <<EOF > ./posteriors
    Id                 Name  Chrom     Position     A1     A2        A1Frq     A1Effect           SE            PIP  LastSampleEff
     1             rs175146     22     17309296      A      G     0.890000    -0.001149     0.005652     0.61250019       0.004332
     2           rs17207051     22     18181984      T      G     0.980000    -0.001620     0.008414     0.65374994      -0.005966
     3            rs5992124     22     18334573      T      C     0.970000     0.003590     0.009958     0.66124958       0.016556
     4             rs361991     22     19263892      T      C     0.520000    -0.000736     0.004496     0.61625034      -0.000928
     5            rs5993853     22     19879938      T      C     0.310000     0.006300     0.007726     0.74750000       0.008830
     6           rs11914240     22     22065251      C      T     0.800000     0.016805     0.010989     0.92249984       0.018251
     7            rs7286558     22     22180183      C      T     0.950000    -0.000941     0.006677     0.63374960       0.013190
     8             rs738865     22     22400423      C      T     0.950000    -0.002468     0.007005     0.64875025      -0.004128
     9            rs5757691     22     22609828      A      G     0.170000    -0.003114     0.006611     0.64875013      -0.002091
EOF

cat <<EOF > ./mapfile_noNA
b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
22:17309296	22:16828406	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:18181984	22:17699218	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
22:18334573	22:17851807	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
22:19263892	22:19276369	rs361991	T	C	rs361991	C	T	rs361991	C	T
22:19879938	22:19892415	rs5993853	T	C	rs5993853	T	C	rs5993853	C	T
22:22065251	22:21710962	rs11914240	C	T	rs11914240	T	C	rs11914240	T	C
22:22180183	22:21825894	rs7286558	C	T	rs7286558	T	C	rs7286558	T	C
22:22400423	22:22046025	rs738865	C	T	rs738865	T	C	rs738865	T	C
22:22609828	22:22255427	rs5757691	A	G	rs5757691	A	G	rs5757691	A	G
EOF

cat <<EOF > ./expected-result_1.tsv
posteriorID	EA	Effect	genoID
rs175146	A	-0.001149	rs175146
rs17207051	T	-0.001620	rs17207051
rs5992124	T	0.003590	rs5992124
rs361991	T	-0.000736	rs361991
rs5993853	T	0.006300	rs5993853
rs11914240	C	0.016805	rs11914240
rs7286558	C	-0.000941	rs7286558
rs738865	C	-0.002468	rs738865
rs5757691	A	-0.003114	rs5757691
EOF

_run_script ${snp_posteriors_cols} ${map_from_to} "true"

#---------------------------------------------------------------------------------
# benchmark
_setup "base_case_benchmark"

snp_posteriors_cols="6,7,9"
map_from_to="3,6"

cat <<EOF > ./posteriors
CHR	POS	CHR	POS	0	RSID	EffectAllele	OtherAllele	B	SE	Z	P	OR	Neff	CaseN	ControlN	EAF	INFO	N
22	17309296	22	16828406	7412256	rs175146	A	G	0.006499	0.0185	0.350318	0.7261	1.00652	127206.851988	54675	76017	0.06	0.925	111487
22	18181984	22	17699218	1576295	rs17207051	T	G	-0.016404	0.0125	-1.307627	0.191	0.98373	150092.970769	64952	88856	0.84	0.809	111487
22	23609726	22	23267539	7357723	rs78919016	A	G	-0.046599	0.0327	-1.426237	0.1538	0.95447	126819.097462	54385	76025	0.97	0.721	111487
EOF

cat <<EOF > ./mapfile_noNA
b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
22:17309296	22:16828406	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:18181984	22:17699218	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
22:18334573	22:17851807	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
EOF

cat <<EOF > ./expected-result_1.tsv
posteriorID	EA	Effect	genoID
rs175146	A	0.006499	rs175146
rs17207051	T	-0.016404	rs17207051
EOF

_run_script ${snp_posteriors_cols} ${map_from_to} "false"

