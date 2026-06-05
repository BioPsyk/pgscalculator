# pgscalculator v2.1 Testing

This folder contains test commands and scenarios for pgscalculator v2.1.

## Quick Start - v2.1 Feature Testing

### Test New Features (INFO/MAF Reference Files)

```bash
cd /home/jesgaaopen/ibp_pipeline_pgscalculator/pgscalculator/test-v2

# Set up all test scenarios
./test_scenarios.sh all

# Run scenario 1 (no reference files - compute MAF from genotypes)
../pgscalculator.sh --config out_test_scenario1/config.yaml --steps prep --sbatch
```

### Test Scenarios

| Scenario | Description | Expected Behavior |
|----------|-------------|-------------------|
| 1 | No INFO/MAF files | MAF computed from genotypes, INFO filter skipped |
| 2 | With INFO file | INFO filter applied, variants < threshold filtered |
| 3 | With MAF file | MAF from provided file, not computed |

### What's New in v2.1

1. **EAF from LD Reference**: `prep-ldref` extracts allele frequencies from LD reference
2. **INFO/MAF Filtering**: `prep-inclusion-list` applies filters with reference files
3. **New Output Format**: `scores.tsv.gz` with ALLELE_CT, N_VARIANTS columns
4. **Augmented Sumstat**: `sumstat_augmented.tsv.gz` with posterior effects
5. **Finalize Step**: `finalize-output` creates final outputs and run summary

### Verifying New Features

After prep completes:

```bash
# Check EAF extraction from LD reference
head out_test_scenario1/prep/references/ldref_eaf.tsv

# Check MAF computation (scenario 1 only)
head out_test_scenario1/prep/references/maf_computed.tsv

# Check inclusion list
wc -l out_test_scenario1/prep/inclusion_list/variant_inclusion_list.tsv
```

After full pipeline:

```bash
# Check new output format
zcat out_test_scenario1/sumstats/*/scores.tsv.gz | head

# Check augmented sumstat
zcat out_test_scenario1/sumstats/*/sumstat_augmented.tsv.gz | head

# Check run summary
cat out_test_scenario1/sumstats/*/details/run_summary.txt
```

---

## Legacy Tests

The original tests from `test-commands-v2.sh` are still available:

| Test | Description | Use Case |
|------|-------------|----------|
| TEST 1 | Normal run (full pipeline) | Verify full workflow works |
| TEST 2 | Reduced sumstat | Quick debugging (~30 variants/chr) |
| TEST 4 | Step-by-step execution | Test modular v2 features |
| TEST 5 | Interactive development | Debugging inside container |
| TEST 6 | Batch processing | Multiple sumstats from list |

## Test Outputs

Test outputs are stored in subdirectories:

| Directory | Description |
|-----------|-------------|
| `out_test_scenario1/` | v2.1 - No reference files |
| `out_test_scenario2/` | v2.1 - With INFO file |
| `out_test_scenario3/` | v2.1 - With MAF file |
| `out_test_normal/` | Legacy - Full pipeline runs |
| `out_test_reduced/` | Legacy - Reduced sumstat runs |
| `out_test_stepwise/` | Legacy - Step-by-step runs |

## Tips

- Use **Scenario 1** to test basic functionality without reference files
- Use **Scenario 2/3** to verify INFO/MAF filtering works
- Add `-d` flag for verbose debug output
- Check `.out` and `.err` files for SLURM job logs

## Config Reference

### Minimal Config (Scenario 1)

```yaml
ld_reference: /path/to/band_ukb_10k_hm3
genotypes: /path/to/genotypes
genotype_manifest: /path/to/manifest.txt
outdir: /path/to/output

chromosomes: 21-22
whichn: totalN
score_columns: 1 2 5

sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  threads: 4
```

### With Reference Files (Scenario 2/3)

```yaml
# Add to above config:
references:
  info_file: /path/to/info_scores.tsv  # Optional
  maf_file: /path/to/maf.tsv           # Optional

filters:
  info_threshold: 0.8
  maf_threshold: 0.01
```

## See Also

- [README.md](../README.md) - Full v2 documentation
- [pipeline-output-spec.md](../docs/pipeline-output-spec.md) - Output format specification
- [pipeline-output-plan.md](../docs/pipeline-output-plan.md) - Implementation details

