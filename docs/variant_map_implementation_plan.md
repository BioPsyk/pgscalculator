% Variant Map Implementation Plan (Dual-Build Support)

## Overview

Implementation changes to support the dual-build variant mapping strategy from
docs/variant_map_mapping_plan.md:

1. Support GRCh38 genotype datasets via pre-augmented LD reference
2. Simplify sumstat processing by matching on native GRCh38 coordinates
3. Pre-compute liftover once (in prep-ldref) and cache for reuse

---

## Phase 1: Config Updates

### Files to modify
- config.template.yaml
- bin/lib/config_parser.sh

### New config options

```yaml
# Genome build of the genotype files (default: GRCh37)
genotype_build: GRCh37  # or GRCh38

# Liftover reference file (pre-sorted on col1)
liftover_reference: /path/to/references/liftover/dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz
```

### Validation logic

```bash
CFG_GENOTYPE_BUILD="${CFG_GENOTYPE_BUILD:-GRCh37}"

# Liftover reference is always required (for dual-position mapfile)
CFG_LIFTOVER_REFERENCE="${CFG_LIFTOVER_REFERENCE:-/path/to/default/liftover_reference.txt.gz}"

if [[ ! -f "${CFG_LIFTOVER_REFERENCE}" ]]; then
    log_error "liftover_reference file not found: ${CFG_LIFTOVER_REFERENCE}"
    exit 1
fi
```

---

## Phase 2: Augmented LD Reference (prep-ldref)

### File: bin/lib/steps/prep_ldref.sh

### Current behavior
- Extracts RSIDs from LD reference .info files
- Outputs chr{N}_ld_rsids with: chr:pos_b37, a1, a2, rsid

### New behavior
- Require liftover reference file (error if missing)
- After extraction, join with liftover reference to add pos_b38
- Output: prep/ldref_augmented/chr{N}_ld_augmented.tsv
- Cache: skip augmentation if augmented files already exist

### New function: augment_ldref_with_liftover()

Logic:
1. ensure_dir augmented_dir
2. Check cache - skip if all chr files present
3. For each chromosome:
   - Join LD ref (pos_b37) with liftover (pos_b37 -> pos_b38)
   - Output: pos_b37, pos_b38, ldref_a1, ldref_a2, ldref_rsid

---

## Phase 3: Updated Prep Mapfile Schema

### File: bin/lib/steps/prep_inclusion_list.sh

### Current schema
```
chr  pos  geno_snpid  geno_a1  geno_a2  ldref_snpid  ldref_a1  ldref_a2  ldref_a2freq
```

### New schema
```
chr  pos_b37  pos_b38  geno_snpid  geno_a1  geno_a2  ldref_snpid  ldref_a1  ldref_a2  ldref_a2freq
```

### Key changes

1. Load augmented LD reference instead of plain chr{N}_ld_rsids
2. Match genotypes based on genotype_build:
   - GRCh37: match on chr:pos_b37 + alleles
   - GRCh38: match on chr:pos_b38 + alleles
3. Both positions come from augmented LD ref (liftover already done)
4. Allele matching includes strand flip support (A/T <-> T/A, C/G <-> G/C)

---

## Phase 4: Simplified Format-Sumstat

### File: bin/lib/steps/format_sumstat.sh

### Current behavior
- Reads both cleaned_GRCh38.gz and cleaned_GRCh37.gz
- Pastes them together to add GRCh37 coordinates

### New behavior
- Read only cleaned_GRCh38.gz
- No paste needed (variant_map has both positions)
- Simple chromosome split

---

## Phase 5: Filter-Variants Matching on pos_b38 + Alleles

### File: bin/lib/steps/filter_variants.sh

### Current behavior
- Matches sumstat to mapfile using chr:pos + alleles (with strand flip support)

### New behavior
- Match sumstat using chr:pos_b38 + alleles (sumstat is GRCh38)
- Mapfile has explicit pos_b37 and pos_b38 columns
- Allele matching with strand flip support unchanged

### Output mapfile header
```
chr  pos_b37  pos_b38  sumstat_snpid  sumstat_effect  sumstat_other  geno_snpid  geno_a1  geno_a2  ldref_snpid  ldref_a1  ldref_a2  ldref_a2freq
```

---

## Phase 6: Liftover Reference Preparation

### Location
```
references/liftover/
  dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz  # sorted on col1 (b37)
  dbsnp_cleansumstat_reference_GRCh38_GRCh37.txt.gz  # sorted on col1 (b38)
```

### Format (space-delimited, LC_ALL=C sorted on column 1)
```
10:1045940 10:1000000 rs1831596373 A C
```
- Column 1: chr:pos of source build (sort key)
- Column 2: chr:pos of target build
- Column 3: rsid
- Column 4-5: alleles

---

## Implementation Order

1. Config changes (no dependencies)
2. prep-ldref augmentation (requires liftover ref files)
3. prep-inclusion-list dual positions (depends on Step 2)
4. format-sumstat simplification (no dependencies)
5. filter-variants pos_b38 matching (depends on Step 3)
6. Testing (all steps)

---

## Backward Compatibility

- Default genotype_build: GRCh37 - no changes for existing workflows
- Liftover reference is required - prep-ldref fails if file is missing
- New mapfile columns added; downstream steps updated to use dual positions

---

## Testing Checklist

- [ ] Config parser accepts new options
- [ ] prep-ldref creates augmented files when liftover available
- [ ] prep-ldref skips augmentation if cached
- [ ] prep-inclusion-list produces mapfile with pos_b37 and pos_b38
- [ ] prep-inclusion-list matches correctly for GRCh37 genotypes
- [ ] prep-inclusion-list matches correctly for GRCh38 genotypes
- [ ] format-sumstat works without GRCh37 paste
- [ ] filter-variants matches sumstat on chr:pos_b38 + alleles
- [ ] Final scores correct for GRCh37 genotype test data
- [ ] Final scores correct for GRCh38 genotype test data
- [ ] SLURM array jobs work with new schema
- [ ] Output variant_map.tsv.gz contains both position columns

---

## Files Changed Summary

| File | Changes |
|------|---------|
| config.template.yaml | Add genotype_build, liftover_reference |
| bin/lib/config_parser.sh | Parse new config options |
| bin/lib/steps/prep_ldref.sh | Add augment_ldref_with_liftover() |
| bin/lib/steps/prep_inclusion_list.sh | Update mapfile schema, matching logic |
| bin/lib/steps/format_sumstat.sh | Remove GRCh37 paste, simplify to split only |
| bin/lib/steps/filter_variants.sh | Match on chr:pos_b38 + alleles, update column indices |
