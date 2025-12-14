# Output Implementation Plan

This document outlines the changes needed to produce the specified output structure from the current v2.1 pipeline.

## Current State vs Target

### Current Output Structure

```
{outdir}/
├── config.yaml
├── prep/
│   ├── genotypes/
│   ├── ldref/
│   ├── inclusion_list/
│   └── variant_map.tsv
└── sumstat_{name}/           # ← Currently named this way
    ├── formatted/
    ├── filtered/
    ├── posteriors/
    ├── posteriors_mapped/
    ├── scores/
    └── scores_combined/
        └── merged.sscore
```

### Target Output Structure

```
{outdir}/
├── config.yaml
├── prep/
│   └── ... (same)
└── sumstats/                 # ← Grouped under sumstats/
    └── {name}/
        ├── scores.tsv.gz           # Final score
        ├── sumstat_augmented.tsv.gz # Augmented sumstat
        ├── details/
        ├── qc/
        └── intermediates/
```

---

## Implementation Tasks

### Phase 1: Directory Structure Changes

**Task 1.1: Change sumstat output path**

| Current | Target |
|---------|--------|
| `{outdir}/sumstat_{name}/` | `{outdir}/sumstats/{name}/` |

Files to modify:
- `pgscalculator-v2.sh` - Update path construction
- `bin/lib/steps/*.sh` - Update all step scripts

**Task 1.2: Reorganize intermediate files**

Move working files under `intermediates/`:

```bash
# Current                          # Target
posteriors/                   →    intermediates/posteriors/
posteriors_mapped/            →    intermediates/posteriors_mapped/
scores/                       →    intermediates/scores/
scores_combined/              →    intermediates/scores_combined/
formatted/                    →    intermediates/formatted/
filtered/                     →    intermediates/filtered/
```

---

### Phase 2: Final Output Generation

**Task 2.1: Generate `scores.tsv.gz`**

Modify `combine_scores.sh` to produce the final format:

```bash
# Current output: merged.sscore
IID     SCORE1_SUM

# Target output: scores.tsv.gz
IID     SCORE_SUM     ALLELE_CT     N_VARIANTS
```

Implementation:
1. Collect ALLELE_CT and variant counts during plink2 scoring
2. Add header with metadata
3. Gzip the output

**Task 2.2: Generate `sumstat_augmented.tsv.gz`**

Create new step: `finalize_output.sh`

This combines:
- Original sumstat (formatted)
- Posterior effects (from posteriors_mapped)
- Variant mapping info (from prep)

```bash
# Join: formatted_sumstat + posteriors + variant_map
awk '
  # Load posteriors
  FNR==NR && ARGIND==1 { post[$1] = $4"\t"$5; next }
  # Load variant map  
  FNR==NR && ARGIND==2 { geno[$1] = $2; next }
  # Process sumstat
  {
    rsid = $3  # RSID column
    print $0, geno[rsid], post[rsid]
  }
' posteriors_combined.tsv variant_map.tsv formatted_sumstat.tsv | gzip > sumstat_augmented.tsv.gz
```

**Task 2.3: Generate `variant_map.tsv.gz`**

Create per-sumstat variant map from prep data + sumstat-specific info:

```bash
# Add sumstat-specific columns (IN_ANALYSIS, etc.)
awk '...' prep/variant_map.tsv intermediates/filtered/* > variant_map.tsv.gz
```

---

### Phase 3: QC Generation

**Task 3.1: Add QC plotting step**

Create new step: `generate_qc.sh`

Uses R or Python to generate:
- `chrall_effect_vs_beta.png` - Scatter: posterior effect vs original beta
- `chrall_position_vs_effect.png` - Manhattan-style: position vs posterior effect
- Per-chromosome versions
- `variant_counts.png` - Bar chart: variants at each step
- `qc_summary.txt` - Numeric stats

Dependencies:
- R with ggplot2, or
- Python with matplotlib

**Task 3.2: QC summary statistics**

> **Note:** The exact level of detail for QC summary is TBD. Options:
> - Per-sumstat summary
> - Per-chromosome breakdown  
> - Combined across all sumstats
> 
> Keeping as placeholder for now.

```bash
# qc_summary.txt (example format - subject to change)
Variants in sumstat:          1,234,567
Variants after INFO filter:   1,100,000
Variants after MAF filter:      900,000
Variants in LD reference:       800,000
Variants in genotypes:          750,000
Variants in final analysis:     700,000

Chromosomes processed:        22
Samples scored:               5,000
Mean score:                   0.0012
SD score:                     0.0045
```

---

### Phase 4: Details/Config Preservation

**Task 4.1: Copy config to details/**

```bash
mkdir -p sumstats/{name}/details/
cp config.yaml sumstats/{name}/details/
echo "Run completed: $(date)" > sumstats/{name}/details/run_log.txt
```

Note: sbayesR parameters are already in `config.yaml`, so no separate file needed.

---

## Step-by-Step Changes

### Step 1: Update `combine_scores.sh`

```diff
- output_file="${scores_dir}/../scores_combined/merged.sscore"
+ output_file="${outdir}/sumstats/${sumstat}/scores.tsv"

- echo -e "IID\tSCORE1_SUM" > "$output_file"
+ echo -e "IID\tSCORE_SUM\tALLELE_CT\tN_VARIANTS" > "$output_file"
```

Also: collect allele counts and variant counts from plink2 output.

### Step 2: Create `finalize_output.sh`

New step that runs after `combine_scores`:

1. Combine all posteriors into single file
2. Join with formatted sumstat
3. Add variant mapping info
4. Generate `sumstat_augmented.tsv.gz`
5. Copy config to details/
6. Generate `variant_map.tsv.gz`

### Step 3: Create `generate_qc.sh`

New optional step:

1. Load posteriors and sumstat
2. Generate QC plots
3. Calculate summary statistics
4. Write `qc_summary.txt`

### Step 4: Update path references

Update all step scripts to use new paths:
- `{sumstat_dir}/posteriors/` → `{sumstat_dir}/intermediates/posteriors/`
- etc.

---

## Priority Order

1. **High Priority** (Core functionality)
   - `scores.tsv.gz` generation (Task 2.1)
   - Directory restructure (Task 1.1, 1.2)
   - `sumstat_augmented.tsv.gz` (Task 2.2)

2. **Medium Priority** (Completeness)
   - Config copy to details/ (Task 4.1)
   - `variant_map.tsv.gz` (Task 2.3)

3. **Low Priority** (Nice to have)
   - QC plots (Task 3.1)
   - QC summary (Task 3.2)

---

## Estimated Effort

| Task | Effort | Files Changed |
|------|--------|---------------|
| Directory restructure | 2 hours | ~10 files |
| scores.tsv.gz | 1 hour | combine_scores.sh |
| sumstat_augmented.tsv.gz | 2 hours | new finalize_output.sh |
| variant_map.tsv.gz | 1 hour | finalize_output.sh |
| Config copy to details/ | 15 min | finalize_output.sh |
| QC generation (TBD) | 3 hours | new generate_qc.sh + R/Python |

**Total: ~9 hours**

---

## Testing Plan

1. Run full pipeline on test sumstat
2. Verify all expected files exist
3. Validate file formats match specification
4. Compare scores with v1 output (should be identical)
5. Check QC plots render correctly

---

## Rollout

1. Implement in feature branch
2. Test on 2-3 sumstats
3. Compare with v1 outputs
4. Merge to main
5. Update documentation



