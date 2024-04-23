#!/usr/bin/env bash

set -euo pipefail

test_script="variant_map_for_sbayesr"
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
  gbuild="$1"
  lbuild="$2"

  "${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "${gbuild}" "${lbuild}" ./observed-results.tsv

  _check_results ./observed-results.tsv ./expected-result_1.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# base case

_setup "base_case gb37 lb37"

cat <<EOF > ./ss2
22:18334573	22:17851807	T	C	rs5992124
22:17309296	22:16828406	A	G	rs175146
22:19263892	22:19276369	T	C	rs361991
22:18181984	22:17699218	T	G	rs17207051
22:17685648	22:17204758	T	C	rs71328205
22:17434931	22:16954041	C	T	rs8138892
22:17697944	22:17217054	G	A	rs78872431
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18181984	G	T	rs17207051
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
22:17434931	T	C	rs8138892
22:17685648	C	T	rs71328205
22:17697944	A	G	rs78872431
EOF

# needs to be sorted
cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs17207051
rs175146
rs361991
rs5992124
rs71328205
rs78872431
EOF

cat <<EOF > ./ld2
22:17309296	G	A	rs175146
22:18181984	G	T	rs17207051
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
EOF

cat <<EOF > ./expected-result_1.tsv
b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
22:17309296	22:16828406	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:17685648	22:17204758	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
22:17697944	22:17217054	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
22:18181984	22:17699218	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
22:18334573	22:17851807	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
22:19263892	22:19276369	rs361991	T	C	rs361991	C	T	rs361991	C	T
EOF

_run_script "37" "37"

#---------------------------------------------------------------------------------
# base case gb38

_setup "base_case gb38 lb38"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
22:19276369	22:19263892	T	C	rs361991
22:17699218	22:18181984	T	G	rs17207051
22:17204758	22:17685648	T	C	rs71328205
22:16954041	22:17434931	C	T	rs8138892
22:17217054	22:17697944	G	A	rs78872431
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18181984	G	T	rs17207051
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
22:17434931	T	C	rs8138892
22:17685648	C	T	rs71328205
22:17697944	A	G	rs78872431
EOF

# needs to be sorted
cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs17207051
rs175146
rs361991
rs5992124
rs71328205
rs78872431
EOF

cat <<EOF > ./ld2
22:17309296	G	A	rs175146
22:18181984	G	T	rs17207051
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
EOF

cat <<EOF > ./expected-result_1.tsv
b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
22:16828406	22:17309296	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:17204758	22:17685648	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
22:17217054	22:17697944	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
22:17699218	22:18181984	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
22:17851807	22:18334573	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
22:19276369	22:19263892	rs361991	T	C	rs361991	C	T	rs361991	C	T
EOF

_run_script "38" "38"

#---------------------------------------------------------------------------------
# base case gb38

_setup "gb38 lb37"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
22:19276369	22:19263892	T	C	rs361991
22:17699218	22:18181984	T	G	rs17207051
22:17204758	22:17685648	T	C	rs71328205
22:16954041	22:17434931	C	T	rs8138892
22:17217054	22:17697944	G	A	rs78872431
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18181984	G	T	rs17207051
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
22:17434931	T	C	rs8138892
22:17685648	C	T	rs71328205
22:17697944	A	G	rs78872431
EOF

# needs to be sorted
cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs17207051
rs175146
rs361991
rs5992124
rs71328205
rs78872431
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17699218	G	T	rs17207051
22:17851807	C	T	rs5992124
22:19276369	C	T	rs361991
EOF

cat <<EOF > ./expected-result_1.tsv
b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
22:16828406	22:17309296	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:17204758	22:17685648	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
22:17217054	22:17697944	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
22:17699218	22:18181984	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
22:17851807	22:18334573	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
22:19276369	22:19263892	rs361991	T	C	rs361991	C	T	rs361991	C	T
EOF

_run_script "38" "37"

#---------------------------------------------------------------------------------
# base case gb38

_setup "gb38 lb37 - NAs in ss"

cat <<EOF > ./ss2
NA:NA	22:17697944	G	A	rs78872431
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
22:19276369	22:19263892	T	C	rs361991
22:17699218	22:18181984	T	G	rs17207051
22:17204758	22:17685648	T	C	rs71328205
22:16954041	22:17434931	C	T	rs8138892
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18181984	G	T	rs17207051
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
22:17434931	T	C	rs8138892
22:17685648	C	T	rs71328205
22:17697944	A	G	rs78872431
EOF

# needs to be sorted
cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs17207051
rs175146
rs361991
rs5992124
rs71328205
rs78872431
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17699218	G	T	rs17207051
22:17851807	C	T	rs5992124
22:19276369	C	T	rs361991
EOF

cat <<EOF > ./expected-result_1.tsv
b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
22:16828406	22:17309296	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:17204758	22:17685648	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
22:17699218	22:18181984	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
22:17851807	22:18334573	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
22:19276369	22:19263892	rs361991	T	C	rs361991	C	T	rs361991	C	T
NA:NA	22:17697944	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
EOF

_run_script "38" "37"

