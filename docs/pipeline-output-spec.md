# Output Specification v2.1

This document defines the expected output structure for pgscalculator v2.1, based on the original v1 outputs with improvements for the modular workflow.

## Overview

The output is organized into:
1. **Prep outputs** - Reusable across all sumstats (run once)
2. **Per-sumstat outputs** - Generated for each trait

## Directory Structure

```
{outdir}/
├── config.yaml                    # Copy of config used
├── prep/                          # Reusable prep outputs
│   ├── genotypes/
│   ├── ldref/
│   ├── inclusion_list/
│   └── variant_map.tsv
│
└── sumstats/
    └── {sumstat_name}/            # Per-sumstat outputs
        ├── scores.tsv.gz          # FINAL: Main PGS scores
        ├── scores_benchmark.tsv.gz # FINAL: Benchmark scores (optional)
        ├── sumstat_augmented.tsv.gz # FINAL: Sumstat with posteriors
        ├── details/               # Config/params used
        ├── qc/                    # QC plots
        └── intermediates/         # Working files (optional cleanup)
```

---

## Final Output Files

### 1. `scores.tsv.gz` - Main PGS Scores

The primary output: polygenic scores for all individuals.

| Column | Type | Description |
|--------|------|-------------|
| IID | string | Sample identifier |
| SCORE_SUM | float | Sum of polygenic effects across all chromosomes |
| ALLELE_CT | int | Total allele count (dosages × 2) |
| N_VARIANTS | int | Number of variants used in scoring |

**Example:**
```
IID	SCORE_SUM	ALLELE_CT	N_VARIANTS
NA20775	0.000227	1234724	617362
HG01873	-0.00243	1234712	617356
```

**Note:** v1 included FID, NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG, FILE_SUM. These are derivable:
- FID: same as IID (or extract from genotype files if needed)
- SCORE1_AVG: SCORE_SUM / N_VARIANTS
- FILE_SUM: always equals number of chromosomes (22)

### 2. `scores_benchmark.tsv.gz` - Benchmark Scores (Optional)

Scores using MAF-filtered + LD-pruned **original effects** (not posteriors).
Used for comparison with the posterior-based scores.

**Method:** LD pruning (`--indep-pairwise 250 50 0.25`) + MAF filter (default 0.05)

| Column | Type | Description |
|--------|------|-------------|
| IID | string | Sample identifier |
| ALLELE_CT | int | Total allele count |
| SCORE1_SUM | float | Sum of pruned effects |

**Note:** Run with `--steps benchmark` after `sumstat` step. Not included in default `--all` flow.

**Current output location:** `{sumstat_dir}/benchmark/benchmark.sscore`

### 3. `sumstat_augmented.tsv.gz` - Augmented Summary Statistics

The input sumstat with posterior effects and variant mapping added.

| Column | Type | Description |
|--------|------|-------------|
| CHR | int | Chromosome |
| POS | int | Position (build 37) |
| RSID | string | rs identifier |
| EA | string | Effect allele |
| OA | string | Other allele |
| EAF | float | Effect allele frequency |
| BETA | float | Original effect size |
| SE | float | Standard error |
| P | float | P-value |
| N | int | Sample size |
| INFO | float | Imputation quality (if available) |
| GENO_ID | string | Matched genotype variant ID |
| POST_EFFECT | float | Posterior effect from sbayesR |
| POST_PIP | float | Posterior inclusion probability |
| IN_ANALYSIS | bool | Whether variant was used in final score |

**Example:**
```
CHR	POS	RSID	EA	OA	EAF	BETA	SE	P	N	GENO_ID	POST_EFFECT	POST_PIP	IN_ANALYSIS
22	16554886	rs78922722	G	T	0.018	-0.0008	0.0072	0.65	100000	rs78922722	-0.000003	0.06	Y
```

### 4. `variant_map.tsv.gz` - Full Variant Mapping

Complete crosswalk between sumstat, genotype, and LD reference variants.

| Column | Type | Description |
|--------|------|-------------|
| RSID | string | rs identifier |
| CHR | int | Chromosome |
| POS_B37 | int | Position in GRCh37 |
| POS_B38 | int | Position in GRCh38 |
| EA | string | Effect allele (from sumstat) |
| OA | string | Other allele (from sumstat) |
| GENO_ID | string | Matched genotype ID |
| GENO_A1 | string | Genotype allele 1 |
| GENO_A2 | string | Genotype allele 2 |
| LD_ID | string | LD reference ID |
| IN_LD | bool | Present in LD reference |
| IN_GENO | bool | Present in genotypes |
| IN_ANALYSIS | bool | Used in analysis |

---

## Supporting Output Files

### 5. `details/` - Run Configuration

```
details/
├── config.yaml          # Full config used
├── sbayesr_params.txt   # sbayesR parameters
└── run_log.txt          # Execution summary
```

### 6. `qc/` - Quality Control Plots

```
qc/
├── chrall_effect_vs_beta.png     # Posterior vs original effect
├── chrall_position_vs_effect.png # Manhattan-style effect plot
├── chr{1-22}_effect_vs_beta.png  # Per-chromosome versions
├── chr{1-22}_position_vs_effect.png
├── variant_counts.png            # Variants per step
└── qc_summary.txt                # Numeric QC stats
```

### 7. `intermediates/` - Working Files

These can be cleaned up after successful completion:

```
intermediates/
├── formatted/
│   ├── sumstat_formatted.tsv.gz
│   └── chr{1-22}_formatted.tsv
├── filtered/
│   ├── sumstat_filtered.tsv.gz
│   └── chr{1-22}_filtered.tsv
├── posteriors/
│   ├── chr{1-22}.snpRes
│   └── work_chr{1-22}/           # sbayesR working files
├── posteriors_mapped/
│   ├── chr{1-22}.snpRes
│   └── rsid_to_genoid.tsv
├── scores/
│   ├── chr{1-22}.sscore
│   └── work_chr{1-22}/
└── posteriors_combined.tsv       # All chromosomes combined
```

---

## Prep Outputs (Reusable)

### `prep/references/` - QC Reference Files

These files enable INFO and MAF filtering during the inclusion list creation.
**Key format:** Uses genotype variant IDs so users can create these from imputation output.

#### `info_scores.tsv` (User-provided)

| Column | Type | Description |
|--------|------|-------------|
| GENO_ID | string | Genotype variant ID (e.g., `chr1:12345:A:G` or `rs12345`) |
| INFO | float | Imputation quality score (0-1) |

```
GENO_ID	INFO
chr1:12345:A:G	0.95
chr1:23456:T:C	0.82
rs789012	0.99
```

**Source:** Extract from imputation server output (e.g., Michigan, TOPMed)

#### `maf.tsv` (Can be computed or user-provided)

| Column | Type | Description |
|--------|------|-------------|
| GENO_ID | string | Genotype variant ID |
| MAF | float | Minor allele frequency (0-0.5) |

```
GENO_ID	MAF
chr1:12345:A:G	0.15
chr1:23456:T:C	0.02
rs789012	0.35
```

**Source:** 
- Computed from genotypes: `plink2 --freq` → extract ID and ALT_FREQS
- Or user-provided from imputation reference panel

#### `ldref_eaf.tsv` (Auto-generated from LD reference)

Extracted from LD reference `.info` files during `prep-ldref` step.

| Column | Type | Description |
|--------|------|-------------|
| RSID | string | rs identifier (from LD reference) |
| A1 | string | Allele 1 |
| A2 | string | Allele 2 |
| A2Freq | float | Frequency of A2 allele (from UKB) |

```
RSID	A1	A2	A2Freq
rs7287144	A	G	0.71675
rs4010558	G	A	0.707
```

**Source:** LD reference `.info` files (e.g., `band_chr22.ldm.sparse.info`)

**Use case:** Fallback EAF for sumstats when EAF column is missing.

#### EAF Priority for Sumstat Processing

1. **EAF** from sumstat (if available and valid)
2. **ldref_eaf** from LD reference (preferred fallback - same population as LD matrix)
3. ~~**EAF_1KG**~~ from sumstat (avoid - different population, less accurate)

#### Config reference

```yaml
references:
  info_file: /path/to/info_scores.tsv    # Optional (user-provided)
  maf_file: /path/to/maf.tsv             # Optional (can compute from genotypes)
  # ldref_eaf is auto-generated during prep-ldref
  
filters:
  info_threshold: 0.8
  maf_threshold: 0.01
```

### `prep/genotypes/`

```
genotypes/
├── snplist_sorted                # All variant IDs, sorted
└── chr{1-22}_pvar_fmt            # Per-chromosome variant lists
```

### `prep/ldref/`

```
ldref/
├── ld_rsids_all                  # All LD reference RSIDs
└── chr{1-22}_ld_rsids            # Per-chromosome LD RSIDs
```

### `prep/inclusion_list/`

```
inclusion_list/
├── variant_inclusion_list.tsv    # Final inclusion list
└── chr{1-22}_variant_map         # Per-chromosome mappings
```

### `prep/variant_map.tsv`

Master variant crosswalk (rsid <-> genotype_id).

---

## Comparison with v1 Output

| v1 File | v2 Equivalent | Changes |
|---------|--------------|---------|
| `main_raw_score_all.gz` | `scores.tsv.gz` | Simplified columns |
| `bench_raw_score_all.gz` | `scores_benchmark.tsv.gz` | Optional |
| `augmented_sumstat.gz` | `sumstat_augmented.tsv.gz` | Added POST_PIP, IN_ANALYSIS |
| `variant_map.gz` | `variant_map.tsv.gz` + `prep/variant_map.tsv` | Split for reuse |
| `details/` | `details/` | Same |
| `extra/raw_posteriors_chrall` | `intermediates/posteriors_combined.tsv` | Combined |
| `extra/raw_maf_chrall` | Removed | Info in augmented_sumstat |
| `qc/` | `qc/` | Same |
| `pipeline_info/` | Removed | Not using Nextflow |
| `intermediates/` | `intermediates/` | Reorganized |

---

## File Sizes (Typical)

| File | Typical Size | Notes |
|------|--------------|-------|
| `scores.tsv.gz` | 10-50 KB | Per sample, small |
| `scores_benchmark.tsv.gz` | 10-50 KB | Same |
| `sumstat_augmented.tsv.gz` | 50-500 MB | All variants |
| `variant_map.tsv.gz` | 10-50 MB | All variants |
| `intermediates/` | 500 MB - 2 GB | Can be cleaned |
| `qc/` | 5-20 MB | PNG plots |

---

## Cleanup Options

After successful completion:

```bash
# Remove working files (keep final outputs)
rm -rf sumstats/{name}/intermediates/

# Keep only essential files
# - scores.tsv.gz
# - sumstat_augmented.tsv.gz
# - variant_map.tsv.gz
# - qc/
```

---

## Validation Checklist

A successful run should have:

- [ ] `scores.tsv.gz` exists and has correct number of samples
- [ ] `sumstat_augmented.tsv.gz` exists with POST_EFFECT column
- [ ] QC plots generated for all chromosomes
- [ ] No empty intermediate files
- [ ] SCORE_SUM values are not all zero or NA



