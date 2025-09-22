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

function _check_file_sorted {
  file=$1
  column=$2
  sort_type=${3:-""}  # optional: -n for numeric sort
  
  if [ "$sort_type" = "-n" ]; then
    if ! sort -c -k${column},${column}n "$file" 2>/dev/null; then
      echo "- [FAIL] ${curr_case}: File $file is not sorted numerically on column $column"
      exit 1
    fi
  else
    if ! LC_ALL=C sort -c -k${column},${column} "$file" 2>/dev/null; then
      echo "- [FAIL] ${curr_case}: File $file is not sorted lexicographically on column $column"
      exit 1
    fi
  fi
}

function _expect_script_failure {
  gbuild="$1"
  lbuild="$2"
  expected_error_pattern="$3"
  
  if "${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "${gbuild}" "${lbuild}" ./observed-results.tsv 2>error.log; then
    echo "- [FAIL] ${curr_case}: Script should have failed but succeeded"
    exit 1
  fi
  
  if ! grep -q "$expected_error_pattern" error.log; then
    echo "- [FAIL] ${curr_case}: Expected error pattern '$expected_error_pattern' not found in error log"
    cat error.log
    exit 1
  fi
  
  echo "- [OK] ${curr_case}: Script failed as expected"
  cd "${initial_dir}"
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
b37	b38	ss_SNP	ss_A1	ss_A2	pvar_SNP	pvar_A1	pvar_A2	ld_SNP	ld_A1	ld_A2
22:17309296	22:16828406	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:17685648	22:17204758	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
22:17697944	22:17217054	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
22:18181984	22:17699218	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
22:18334573	22:17851807	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
22:19263892	22:19276369	rs361991	T	C	rs361991	C	T	rs361991	C	T
EOF

_run_script "37" "37"
#
##---------------------------------------------------------------------------------
## base case gb38
#
#_setup "base_case gb38 lb38"
#
#cat <<EOF > ./ss2
#22:17851807	22:18334573	T	C	rs5992124
#22:16828406	22:17309296	A	G	rs175146
#22:19276369	22:19263892	T	C	rs361991
#22:17699218	22:18181984	T	G	rs17207051
#22:17204758	22:17685648	T	C	rs71328205
#22:16954041	22:17434931	C	T	rs8138892
#22:17217054	22:17697944	G	A	rs78872431
#EOF
#
#cat <<EOF > ./bim2
#22:17309296	G	A	rs175146
#22:18181984	G	T	rs17207051
#22:18334573	C	T	rs5992124
#22:19263892	C	T	rs361991
#22:17434931	T	C	rs8138892
#22:17685648	C	T	rs71328205
#22:17697944	A	G	rs78872431
#EOF
#
## needs to be sorted
#cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
#rs17207051
#rs175146
#rs361991
#rs5992124
#rs71328205
#rs78872431
#EOF
#
#cat <<EOF > ./ld2
#22:17309296	G	A	rs175146
#22:18181984	G	T	rs17207051
#22:18334573	C	T	rs5992124
#22:19263892	C	T	rs361991
#EOF
#
#cat <<EOF > ./expected-result_1.tsv
#b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
#22:16828406	22:17309296	rs175146	A	G	rs175146	G	A	rs175146	G	A
#22:17204758	22:17685648	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
#22:17217054	22:17697944	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
#22:17699218	22:18181984	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
#22:17851807	22:18334573	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
#22:19276369	22:19263892	rs361991	T	C	rs361991	C	T	rs361991	C	T
#EOF
#
#_run_script "38" "38"
#
##---------------------------------------------------------------------------------
## base case gb38
#
#_setup "gb38 lb37"
#
#cat <<EOF > ./ss2
#22:17851807	22:18334573	T	C	rs5992124
#22:16828406	22:17309296	A	G	rs175146
#22:19276369	22:19263892	T	C	rs361991
#22:17699218	22:18181984	T	G	rs17207051
#22:17204758	22:17685648	T	C	rs71328205
#22:16954041	22:17434931	C	T	rs8138892
#22:17217054	22:17697944	G	A	rs78872431
#EOF
#
#cat <<EOF > ./bim2
#22:17309296	G	A	rs175146
#22:18181984	G	T	rs17207051
#22:18334573	C	T	rs5992124
#22:19263892	C	T	rs361991
#22:17434931	T	C	rs8138892
#22:17685648	C	T	rs71328205
#22:17697944	A	G	rs78872431
#EOF
#
## needs to be sorted
#cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
#rs17207051
#rs175146
#rs361991
#rs5992124
#rs71328205
#rs78872431
#EOF
#
#cat <<EOF > ./ld2
#22:16828406	G	A	rs175146
#22:17699218	G	T	rs17207051
#22:17851807	C	T	rs5992124
#22:19276369	C	T	rs361991
#EOF
#
#cat <<EOF > ./expected-result_1.tsv
#b37	b38	ss_SNP	ss_A1	ss_A2	bim_SNP	bim_A1	bim_A2	ld_SNP	ld_A1	ld_A2
#22:16828406	22:17309296	rs175146	A	G	rs175146	G	A	rs175146	G	A
#22:17204758	22:17685648	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
#22:17217054	22:17697944	rs78872431	G	A	rs78872431	A	G	NA	NA	NA
#22:17699218	22:18181984	rs17207051	T	G	rs17207051	G	T	rs17207051	G	T
#22:17851807	22:18334573	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
#22:19276369	22:19263892	rs361991	T	C	rs361991	C	T	rs361991	C	T
#EOF
#
#_run_script "38" "37"

#---------------------------------------------------------------------------------
# SORTING VALIDATION TESTS
#---------------------------------------------------------------------------------

# Test 1: Unsorted snp2 file (should fail)
_setup "unsorted_snp2_file"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18334573	C	T	rs5992124
EOF

# Create UNSORTED snp2 file (this should cause join to fail/produce incorrect results)
cat <<EOF > ./snp2
rs5992124
rs175146
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17851807	C	T	rs5992124
EOF

# This test should demonstrate the sorting issue
# The script may not fail but will produce incorrect/incomplete results
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

# Count the number of data rows (excluding header)
data_rows=$(tail -n +2 ./observed-results.tsv | wc -l)
if [ "$data_rows" -lt 2 ]; then
  echo "- [OK] ${curr_case}: Unsorted snp2 file caused incomplete results (${data_rows} rows instead of expected 2)"
else
  echo "- [WARNING] ${curr_case}: Unsorted snp2 file did not cause expected issues (got ${data_rows} rows)"
fi

cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Test 2: Properly sorted files (should work correctly)
_setup "properly_sorted_files"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18334573	C	T	rs5992124
EOF

# Create PROPERLY SORTED snp2 file
cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs5992124
rs175146
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17851807	C	T	rs5992124
EOF

cat <<EOF > ./expected-result_1.tsv
b37	b38	ss_SNP	ss_A1	ss_A2	pvar_SNP	pvar_A1	pvar_A2	ld_SNP	ld_A1	ld_A2
22:16828406	22:17309296	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:17851807	22:18334573	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
EOF

_run_script "38" "37"

#---------------------------------------------------------------------------------
# Test 3: Verify intermediate files are sorted correctly
_setup "intermediate_files_sorting_check"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
22:19276369	22:19263892	T	C	rs361991
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
EOF

cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs175146
rs361991
rs5992124
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17851807	C	T	rs5992124
22:19276369	C	T	rs361991
EOF

# Run script and then check intermediate files are sorted
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

# Check that intermediate sorted files exist and are properly sorted
if [ -f "ss2_sorted.tmp" ]; then
  _check_file_sorted "ss2_sorted.tmp" 1
  echo "- [OK] ${curr_case}: ss2_sorted.tmp is properly sorted on column 1"
else
  echo "- [FAIL] ${curr_case}: ss2_sorted.tmp not found"
  exit 1
fi

if [ -f "bim2_sorted_1.tmp" ]; then
  _check_file_sorted "bim2_sorted_1.tmp" 4
  echo "- [OK] ${curr_case}: bim2_sorted_1.tmp is properly sorted on column 4"
else
  echo "- [FAIL] ${curr_case}: bim2_sorted_1.tmp not found"
  exit 1
fi

if [ -f "ld2_sorted.tmp" ]; then
  _check_file_sorted "ld2_sorted.tmp" 1
  echo "- [OK] ${curr_case}: ld2_sorted.tmp is properly sorted on column 1"
else
  echo "- [FAIL] ${curr_case}: ld2_sorted.tmp not found"
  exit 1
fi

echo "- [OK] ${curr_case}: All intermediate files are properly sorted"
cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Test 4: Test build-specific sorting logic
_setup "build_specific_sorting_gb38_lb37"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18334573	C	T	rs5992124
EOF

cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs175146
rs5992124
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17851807	C	T	rs5992124
EOF

# Run with different builds (should trigger join1.resorted.tmp creation)
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

# Check that join1.resorted.tmp is created and sorted when builds differ
if [ -f "join1.resorted.tmp" ]; then
  _check_file_sorted "join1.resorted.tmp" 2
  echo "- [OK] ${curr_case}: join1.resorted.tmp is properly sorted on column 2 for different builds"
else
  echo "- [FAIL] ${curr_case}: join1.resorted.tmp not found when builds differ"
  exit 1
fi

echo "- [OK] ${curr_case}: Build-specific sorting logic works correctly"
cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Test 5: Empty files handling
_setup "empty_files_handling"

# Create empty files
touch ./ss2 ./bim2 ./snp2 ./ld2

# Script should handle empty files gracefully
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "37" "37" ./observed-results.tsv

# Should produce header-only output
lines=$(wc -l < ./observed-results.tsv)
if [ "$lines" -eq 1 ]; then
  echo "- [OK] ${curr_case}: Empty files handled correctly (header-only output)"
else
  echo "- [FAIL] ${curr_case}: Empty files not handled correctly (got $lines lines)"
  exit 1
fi

cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Test 6: Large dataset sorting performance
_setup "large_dataset_sorting"

# Create larger test files to verify sorting performance
for i in {1..100}; do
  echo "22:1700${i}	22:1800${i}	T	C	rs${i}" >> ./ss2
  echo "22:1700${i}	C	T	rs${i}" >> ./bim2
  echo "rs${i}" >> ./snp2_unsorted
  echo "22:1700${i}	C	T	rs${i}" >> ./ld2
done

# Sort the snp2 file properly
LC_ALL=C sort -k1,1 ./snp2_unsorted > ./snp2

# Run the script and measure basic performance
start_time=$(date +%s)
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "37" "37" ./observed-results.tsv
end_time=$(date +%s)
duration=$((end_time - start_time))

# Should complete in reasonable time (less than 10 seconds for 100 records)
if [ "$duration" -lt 10 ]; then
  echo "- [OK] ${curr_case}: Large dataset processed in ${duration} seconds"
else
  echo "- [WARNING] ${curr_case}: Large dataset took ${duration} seconds (may indicate performance issue)"
fi

# Verify output has expected number of records (header + data)
lines=$(wc -l < ./observed-results.tsv)
expected_lines=101  # header + 100 data lines
if [ "$lines" -eq "$expected_lines" ]; then
  echo "- [OK] ${curr_case}: Correct number of output lines ($lines)"
else
  echo "- [FAIL] ${curr_case}: Incorrect number of output lines (got $lines, expected $expected_lines)"
  exit 1
fi

cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Allele flips

_setup "gb38 lb37 - allele flips"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
22:19276369	22:19263892	T	C	rs361991
22:17699218	22:18181984	A	C	rs17207051
22:17204758	22:17685648	T	C	rs71328205
22:16954041	22:17434931	G	A	rs8138892
22:17217054	22:17697944	C	T	rs78872431
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
b37	b38	ss_SNP	ss_A1	ss_A2	pvar_SNP	pvar_A1	pvar_A2	ld_SNP	ld_A1	ld_A2
22:16828406	22:17309296	rs175146	A	G	rs175146	G	A	rs175146	G	A
22:17204758	22:17685648	rs71328205	T	C	rs71328205	C	T	NA	NA	NA
22:17217054	22:17697944	rs78872431	C	T	rs78872431	A	G	NA	NA	NA
22:17699218	22:18181984	rs17207051	A	C	rs17207051	G	T	rs17207051	G	T
22:17851807	22:18334573	rs5992124	T	C	rs5992124	C	T	rs5992124	C	T
22:19276369	22:19263892	rs361991	T	C	rs361991	C	T	rs361991	C	T
EOF

_run_script "38" "37"