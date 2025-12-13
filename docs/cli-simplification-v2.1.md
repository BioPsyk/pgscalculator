# CLI Simplification - pgscalculator v2.1

## What Changed (v2.0 → v2.1)

### Removed Options
| Option | Reason |
|--------|--------|
| `-l` | Moved to config: `ld_reference` |
| `-g` | Moved to config: `genotypes` |
| `-f` | Moved to config: `genotype_manifest` |
| `-c` | Replaced by `--config` |
| `-1` | Deprecated, use `--steps` |
| `-2` | Deprecated, use `--steps` |
| `-b` | Removed (auto-handled) |
| `-w` | Removed (auto-handled) |
| `-j` | Removed (singularity default) |

### Kept Options
| Option | Purpose |
|--------|---------|
| `-i` | Sumstat input path (absolute) |
| `-o` | Output directory (optional, overrides config) |
| `--config` | **NEW: Required** config.yaml path |
| `--steps` | Step groups to run |
| `--skip-prep` | Skip prep if done |
| `--chr` | Chromosome range |
| `-d` | Verbose mode |

---

## Current Design (v2.1 - IMPLEMENTED)

### Core Philosophy

1. **Config-first**: Paths go in `config.yaml`, not CLI
2. **Prep once, reuse**: Prep outputs are shared across sumstats
3. **Simple CLI**: Minimal required arguments
4. **sbatch-friendly**: Works naturally with job arrays

### Simplified CLI

```bash
# Run prep steps (once per project)
./pgscalculator-v2.sh --config config.yaml --steps prep

# Run all steps for a sumstat
./pgscalculator-v2.sh --config config.yaml -i /path/to/sumstat_TRAIT

# Run per-sumstat steps (after prep is done)
./pgscalculator-v2.sh --config config.yaml -i /path/to/sumstat_TRAIT --skip-prep

# Run specific steps only
./pgscalculator-v2.sh --config config.yaml -i /path/to/sumstat_TRAIT --steps posteriors,score --skip-prep
```

### Config-Based Approach

All reference paths in `config.yaml`:

```yaml
# config.yaml - Project configuration

# Output directory
outdir: /path/to/project_output

# Reference data paths (required)
ld_reference: /path/to/band_ukb_10k_hm3
genotypes: /path/to/genotypes
genotype_manifest: /path/to/manifest.txt

# Filtering thresholds
info_threshold: 0.8
maf_threshold: 0.01
  
# sbayesR parameters
sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  threads: 6
  seed: 80851
  exclude_mhc: true

# Optional: chromosome subset for testing
# chromosomes: "21-22"
```

### CLI Structure (v2.1)

```bash
./pgscalculator-v2.sh --config <file> [options]

# Required
--config FILE          # Path to config.yaml (required)

# Optional
-i PATH                # Sumstat input path (absolute, required for non-prep)
-o PATH                # Output directory (overrides config)
--steps GROUPS         # Step groups: prep,sumstat,posteriors,score
--skip-prep            # Skip prep if already done
--chr RANGE            # Chromosome range (e.g., "21-22")
-d                     # Verbose/debug output

# Help
-h, --help             # Show help
-v                     # Show version
```

---

## Recommended Workflows

### Workflow 1: Full Project Setup

```bash
# 1. Create project
mkdir my_pgs_project && cd my_pgs_project

# 2. Initialize config
pgscalculator init \
  --genotypes /data/genotypes \
  --ld-reference /data/ld/band_ukb_10k_hm3 \
  --output ./output

# 3. Run prep (submit as sbatch job)
sbatch --mem=10g --time=1:00:00 --wrap="
  pgscalculator prep --config config.yaml
"

# 4. Process sumstats (job array or loop)
for id in 814 815 816; do
  sbatch --mem=20g --time=2:00:00 --wrap="
    pgscalculator run \
      --sumstat /data/sumstats/sumstat_${id} \
      --config config.yaml
  "
done
```

### Workflow 2: Batch Processing with SLURM Array

```bash
# sumstat_list.txt contains one sumstat path per line
# /data/sumstats/sumstat_814
# /data/sumstats/sumstat_815
# ...

sbatch --array=1-100 --mem=20g --time=2:00:00 <<'EOF'
#!/bin/bash
SUMSTAT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sumstat_list.txt)
pgscalculator run --sumstat "$SUMSTAT" --config config.yaml
EOF
```

### Workflow 3: Step-by-Step for Debugging

```bash
# Prep (run once)
pgscalculator prep --config config.yaml

# Format sumstat
pgscalculator run \
  --sumstat /data/sumstats/sumstat_814 \
  --config config.yaml \
  --steps sumstat

# Calculate posteriors only
pgscalculator run \
  --sumstat /data/sumstats/sumstat_814 \
  --config config.yaml \
  --steps posteriors \
  --skip-prep

# Calculate scores only  
pgscalculator run \
  --sumstat /data/sumstats/sumstat_814 \
  --config config.yaml \
  --steps score \
  --skip-prep
```

---

## Implementation Plan

### Phase 1: Clean up v2.0 (Quick Wins)

1. Remove `-1` and `-2` options from help and code
2. Update documentation to emphasize `--steps`
3. Add deprecation warning if `-1`/`-2` used

### Phase 2: Config Migration

1. Create `pgscalculator init` command
2. Move genotype/LD paths to config
3. Keep `-i` (sumstat) and `-o` (output) as CLI for backwards compatibility
4. Update wrapper to read paths from config

### Phase 3: Simplify Container Handling

1. Auto-detect singularity vs docker
2. Build mount paths from config
3. Better error messages for missing mounts

### Phase 4: Enhanced Workflow Support

1. Add `pgscalculator status` improvements
2. Job dependency tracking
3. Resume failed runs automatically

---

## Output Directory Structure (v2.1)

```
project_output/
├── config.yaml              # Project configuration
├── prep/                    # Shared prep outputs (run once)
│   ├── genotypes/
│   ├── ldref/
│   ├── inclusion_list/
│   └── variant_map.tsv
├── sumstat_814/             # Per-sumstat outputs
│   ├── formatted/
│   ├── filtered/
│   ├── posteriors/
│   ├── posteriors_mapped/
│   ├── scores/
│   └── scores_combined/
├── sumstat_815/
│   └── ...
└── logs/
    ├── prep.log
    ├── sumstat_814.log
    └── ...
```

---

## Questions for Discussion

1. **Should `--sumstat` accept a sumstat ID instead of full path?**
   - Could look up in a configured sumstat library
   - Example: `--sumstat 814` → looks in `/data/sumstats/sumstat_814`

2. **Should we support `--from-step` for resuming?**
   - Example: `--from-step calc-posteriors` to skip format/filter

3. **Should container image be configurable per-run or locked to project?**
   - Reproducibility vs flexibility tradeoff

4. **Multi-method support (sbayesR + PRS-CS)?**
   - Could add `--method sbayesr|prscs` option
   - Or separate commands: `pgscalculator sbayesr`, `pgscalculator prscs`

