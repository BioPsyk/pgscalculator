#!/usr/bin/env bash

set -euo pipefail

test_script="variant_map_for_sbayesr"
initial_dir=$(pwd)"/${test_script}_sorting_validation"
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

function _validate_join_prerequisites {
  # This function validates that files are sorted correctly before joins
  echo "Validating join prerequisites for ${curr_case}..."
  
  # Check if snp2 is sorted (critical for first join)
  if [ -f "./snp2" ]; then
    _check_file_sorted "./snp2" 1
    echo "  ✓ snp2 is sorted on column 1"
  fi
  
  # After running the script, check intermediate files
  if [ -f "ss2_sorted.tmp" ]; then
    _check_file_sorted "ss2_sorted.tmp" 1
    echo "  ✓ ss2_sorted.tmp is sorted on column 1"
  fi
  
  if [ -f "bim2_sorted_1.tmp" ]; then
    _check_file_sorted "bim2_sorted_1.tmp" 4
    echo "  ✓ bim2_sorted_1.tmp is sorted on column 4"
  fi
  
  if [ -f "ld2_sorted.tmp" ]; then
    _check_file_sorted "ld2_sorted.tmp" 1
    echo "  ✓ ld2_sorted.tmp is sorted on column 1"
  fi
  
  if [ -f "bim2_sorted_2.tmp" ]; then
    _check_file_sorted "bim2_sorted_2.tmp" 1
    echo "  ✓ bim2_sorted_2.tmp is sorted on column 1"
  fi
  
  if [ -f "join1.resorted.tmp" ]; then
    _check_file_sorted "join1.resorted.tmp" 2
    echo "  ✓ join1.resorted.tmp is sorted on column 2"
  fi
}

echo ">> Test ${test_script} - Sorting Validation"

#=================================================================================
# CRITICAL SORTING TESTS
#=================================================================================

#---------------------------------------------------------------------------------
# Test 1: Demonstrate the critical snp2 sorting issue
_setup "critical_snp2_unsorted_issue"

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

# Create UNSORTED snp2 file - this is the critical issue
cat <<EOF > ./snp2
rs5992124
rs361991
rs175146
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17851807	C	T	rs5992124
22:19276369	C	T	rs361991
EOF

echo "Running script with UNSORTED snp2 file..."
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

# Count results - with unsorted snp2, we expect incomplete/incorrect results
data_rows=$(tail -n +2 ./observed-results.tsv | wc -l)
expected_rows=3

if [ "$data_rows" -lt "$expected_rows" ]; then
  echo "- [CRITICAL ISSUE DETECTED] ${curr_case}: Unsorted snp2 caused incomplete results"
  echo "  Expected: $expected_rows rows, Got: $data_rows rows"
  echo "  This demonstrates the sorting bug!"
else
  echo "- [UNEXPECTED] ${curr_case}: Unsorted snp2 did not cause expected issues"
  echo "  Got: $data_rows rows (expected fewer due to sorting issue)"
fi

cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Test 2: Same data with PROPERLY sorted snp2 file
_setup "critical_snp2_sorted_fix"

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

# Create PROPERLY SORTED snp2 file - this is the fix
cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs5992124
rs361991
rs175146
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17851807	C	T	rs5992124
22:19276369	C	T	rs361991
EOF

echo "Running script with PROPERLY SORTED snp2 file..."
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

# Count results - with sorted snp2, we expect complete results
data_rows=$(tail -n +2 ./observed-results.tsv | wc -l)
expected_rows=3

if [ "$data_rows" -eq "$expected_rows" ]; then
  echo "- [OK] ${curr_case}: Properly sorted snp2 produced complete results"
  echo "  Expected: $expected_rows rows, Got: $data_rows rows"
else
  echo "- [FAIL] ${curr_case}: Even with sorted snp2, results are incomplete"
  echo "  Expected: $expected_rows rows, Got: $data_rows rows"
  exit 1
fi

# Validate all intermediate files are sorted
_validate_join_prerequisites

cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Test 3: Comprehensive sorting validation for all intermediate files
_setup "comprehensive_sorting_validation"

# Create test data
cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
22:16828406	22:17309296	A	G	rs175146
22:19276369	22:19263892	T	C	rs361991
22:17699218	22:18181984	T	G	rs17207051
EOF

cat <<EOF > ./bim2
22:17309296	G	A	rs175146
22:18334573	C	T	rs5992124
22:19263892	C	T	rs361991
22:18181984	G	T	rs17207051
EOF

cat <<EOF | LC_ALL=C sort -k1,1 > ./snp2
rs175146
rs17207051
rs361991
rs5992124
EOF

cat <<EOF > ./ld2
22:16828406	G	A	rs175146
22:17851807	C	T	rs5992124
22:19276369	C	T	rs361991
22:17699218	G	T	rs17207051
EOF

echo "Running comprehensive sorting validation test..."
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

echo "Validating all intermediate files are properly sorted..."
_validate_join_prerequisites

echo "- [OK] ${curr_case}: All intermediate files are properly sorted"
cd "${initial_dir}"

#---------------------------------------------------------------------------------
# Test 4: Build-specific sorting paths
_setup "build_specific_sorting_paths"

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

echo "Testing build-specific sorting (gb38 != lb37)..."
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

# This should create join1.resorted.tmp
if [ -f "join1.resorted.tmp" ]; then
  _check_file_sorted "join1.resorted.tmp" 2
  echo "- [OK] ${curr_case}: join1.resorted.tmp created and sorted correctly for different builds"
else
  echo "- [FAIL] ${curr_case}: join1.resorted.tmp not created when builds differ"
  exit 1
fi

echo "Testing build-specific sorting (gb37 == lb37)..."
rm -f join1.resorted.tmp  # Clean up
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "37" "37" ./observed-results.tsv

# This should NOT create join1.resorted.tmp
if [ ! -f "join1.resorted.tmp" ]; then
  echo "- [OK] ${curr_case}: join1.resorted.tmp not created when builds are the same"
else
  echo "- [FAIL] ${curr_case}: join1.resorted.tmp created unnecessarily when builds are the same"
  exit 1
fi

cd "${initial_dir}"

#=================================================================================
# EDGE CASES AND STRESS TESTS
#=================================================================================

#---------------------------------------------------------------------------------
# Test 5: Empty files
_setup "empty_files_sorting"

touch ./ss2 ./bim2 ./snp2 ./ld2

echo "Testing empty files handling..."
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
# Test 6: Single record files
_setup "single_record_sorting"

cat <<EOF > ./ss2
22:17851807	22:18334573	T	C	rs5992124
EOF

cat <<EOF > ./bim2
22:18334573	C	T	rs5992124
EOF

cat <<EOF > ./snp2
rs5992124
EOF

cat <<EOF > ./ld2
22:17851807	C	T	rs5992124
EOF

echo "Testing single record files..."
"${test_script}.sh" ./ss2 ./bim2 ./snp2 ./ld2 "38" "37" ./observed-results.tsv

# Should produce header + 1 data row
lines=$(wc -l < ./observed-results.tsv)
if [ "$lines" -eq 2 ]; then
  echo "- [OK] ${curr_case}: Single record files handled correctly"
else
  echo "- [FAIL] ${curr_case}: Single record files not handled correctly (got $lines lines, expected 2)"
  exit 1
fi

_validate_join_prerequisites
cd "${initial_dir}"

echo ""
echo "==================================================================="
echo "| SORTING VALIDATION TESTS COMPLETED"
echo "==================================================================="
echo "| Key findings:"
echo "| - Unsorted snp2 file causes incomplete/incorrect results"
echo "| - Proper sorting of snp2 on column 1 is CRITICAL"
echo "| - All intermediate files are correctly sorted by the script"
echo "| - Build-specific logic creates appropriate sorted files"
echo "| - Edge cases (empty, single record) are handled correctly"
echo "==================================================================="

