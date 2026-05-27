## pgscalculator v2.2 smoke tests (LDpred2)

Smoke configs for **sBayesR-only**, **LDpred2-only**, and **both methods**, alongside the v2.1 layout.

Reuses the same **genotypes** and **sBayesR LD reference** paths as `tests/smoke/v2.1-2026-01-11/`.

### Paths to override

| Variable | Purpose |
|----------|---------|
| `SUMSTAT_DIR` | Cleansumstats output folder (one trait subfolder is enough for a quick run) |
| `LDPRED2_LD_DIR` | Directory with `LD_with_blocks_chr*.rds` and `map_hm3_plus.rds` (see `docs/references.md`) |

Default sumstat folder matches v2.1 (`sumstat_5759` on faststorage).

### Run sBayesR only (Day-1 incremental workflow)

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator
export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"
bash tests/smoke/v2.2-2026-05-26/run_local.sh config.sbayesr.yaml
```

### Run LDpred2 only (chr 22 subset, ~2 min interactive with HM3)

```bash
export LDPRED2_LD_DIR="/path/to/references/ld-ldpred2/hm3_plus"
bash tests/smoke/v2.2-2026-05-26/run_local.sh config.ldpred2.yaml
```

### Run both methods in one submission

```bash
bash tests/smoke/v2.2-2026-05-26/run_local.sh config.both.yaml
```

### Incremental run (§6.5)

```bash
# Day 1
bash tests/smoke/v2.2-2026-05-26/run_local.sh config.sbayesr.yaml

# Day 2 — same outdir, add LDpred2 only
export LDPRED2_LD_DIR="/path/to/references/ld-ldpred2/hm3_plus"
bash tests/smoke/v2.2-2026-05-26/run_incremental_day2.sh
```

After both days, check per-sumstat outputs:

- `scores_sbayesr.gz` from Day 1 unchanged
- `scores_ldpred2.gz` from Day 2
- `augmented_sumstat.gz` header includes `postEffect_sbayesr` and `postEffect_ldpred2`

### Method correlation (when both score files exist)

```bash
Rscript tests/smoke/v2.2-2026-05-26/compare_methods.R \
  /path/to/out/sumstats/TRAIT_NAME/scores_sbayesr.gz \
  /path/to/out/sumstats/TRAIT_NAME/scores_ldpred2.gz
```

Expect sample-level Pearson r > 0.6 for well-powered traits.

### SLURM

```bash
bash tests/smoke/v2.2-2026-05-26/submit_slurm.sh config.both.yaml
```
