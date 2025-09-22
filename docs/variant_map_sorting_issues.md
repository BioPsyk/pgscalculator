# Variant Mapping Sorting Issues Analysis

## Overview

The `variant_map_for_sbayesr.sh` script performs multiple join operations on genomic data files. For these joins to work correctly, all input files must be properly sorted on the columns being joined. This document identifies critical sorting issues that can cause join failures or incorrect results.

## Critical Issues Identified

### 1. **Missing Sort for `snp2_sorted` File (CRITICAL)**

**Location**: Line 33 - First join operation
```bash
LC_ALL=C join -1 1 -2 4 -o 2.1 2.2 2.3 2.4 "${snp2_sorted}" bim2_sorted_1.tmp > join0.tmp
```

**Problem**: 
- The `snp2_sorted` file is used directly without being sorted
- The join expects column 1 of `snp2_sorted` to be sorted to match with column 4 of `bim2_sorted_1.tmp`
- This will cause the join to fail or produce incorrect/incomplete results

**Impact**: HIGH - Join will not work correctly, leading to missing or incorrect variant mappings

### 2. **Sorting Operations Analysis**

| Operation | File | Sort Key | Status | Notes |
|-----------|------|----------|--------|-------|
| Line 18/21 | `ss2_sorted.tmp` | Column 1 | ✅ CORRECT | Properly sorted |
| Line 29 | `bim2_sorted_1.tmp` | Column 4 | ✅ CORRECT | Properly sorted |
| Line 30 | `ld2_sorted.tmp` | Column 1 | ✅ CORRECT | Properly sorted |
| Line 34 | `bim2_sorted_2.tmp` | Column 1 | ✅ CORRECT | Re-sorted from join0 output |
| Line 43 | `join1.resorted.tmp` | Column 2 | ✅ CORRECT | Conditional re-sort |
| **MISSING** | `snp2_sorted` | Column 1 | ❌ **MISSING** | **Critical issue** |

### 3. **Join Operations Analysis**

| Join | Line | File 1 | Join Key 1 | File 2 | Join Key 2 | Status |
|------|------|--------|------------|--------|------------|--------|
| Join 0 | 33 | `snp2_sorted` | Col 1 (❌ unsorted) | `bim2_sorted_1.tmp` | Col 4 (✅ sorted) | **BROKEN** |
| Join 1 | 36 | `ss2_sorted.tmp` | Col 1 (✅ sorted) | `bim2_sorted_2.tmp` | Col 1 (✅ sorted) | ✅ OK |
| Join 2a | 41 | `join1.tmp` | Col 1 (✅ sorted) | `ld2_sorted.tmp` | Col 1 (✅ sorted) | ✅ OK |
| Join 2b | 44 | `join1.resorted.tmp` | Col 2 (✅ sorted) | `ld2_sorted.tmp` | Col 1 (✅ sorted) | ✅ OK |

## Required Fixes

### Immediate Fix Required

Add sorting for `snp2_sorted` before the first join:

```bash
# Add after line 14:
LC_ALL=C sort -k1,1 "${snp2_sorted}" > snp2_sorted_fixed.tmp

# Update line 33 to use the sorted file:
LC_ALL=C join -1 1 -2 4 -o 2.1 2.2 2.3 2.4 "snp2_sorted_fixed.tmp" bim2_sorted_1.tmp > join0.tmp
```

## Testing Requirements

To prevent these issues from occurring again, comprehensive unit tests must be implemented to validate:

### 1. **Pre-Join Sorting Validation Tests**
- Verify all input files are properly sorted before join operations
- Test with unsorted input files to ensure they fail gracefully
- Validate sort order matches join requirements

### 2. **Join Operation Integrity Tests**
- Test each join operation with known sorted and unsorted inputs
- Verify join results are complete and correct
- Test edge cases with empty files, single records, and duplicate keys

### 3. **Build-Specific Logic Tests**
- Test different build combinations (b37/b38) that trigger different join paths
- Verify conditional sorting logic works correctly
- Test column swapping logic for build 38

### 4. **Data Integrity Tests**
- Verify allele matching logic works with properly sorted data
- Test complement allele flipping with sorted inputs
- Validate output format and completeness

### 5. **Performance and Memory Tests**
- Test with large files to ensure sorting doesn't cause memory issues
- Verify LC_ALL=C locale setting is consistent across all operations
- Test temporary file cleanup

## Test Data Requirements

Each test should include:
- **Unsorted input files** to verify the script fails appropriately
- **Properly sorted input files** to verify correct operation
- **Mixed scenarios** with some sorted and some unsorted files
- **Edge cases** like empty files, single records, duplicate keys
- **Different build combinations** (b37, b38, mixed)
- **Allele flip scenarios** with proper sorting

## Expected Test Outcomes

1. **Sorting validation tests** should catch the missing sort issue immediately
2. **Join integrity tests** should verify all joins produce expected results
3. **Regression tests** should prevent future sorting issues
4. **Performance tests** should ensure the fixes don't impact performance

## Implementation Priority

1. **HIGH**: Fix the critical missing sort for `snp2_sorted`
2. **HIGH**: Implement sorting validation tests
3. **MEDIUM**: Add comprehensive join integrity tests
4. **MEDIUM**: Add build-specific logic tests
5. **LOW**: Add performance and edge case tests

## Notes

- All sort operations should use `LC_ALL=C` for consistent locale-independent sorting
- Temporary files should be properly cleaned up
- Error handling should be added for sort failures
- Consider adding verbose logging for debugging join operations

