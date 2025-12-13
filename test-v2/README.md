# pgscalculator v2 Testing

This folder contains test commands for pgscalculator v2.

## Quick Start

1. **Edit paths** in `test-commands-v2.sh` to match your environment
2. **Request an interactive node** (for development):
   ```bash
   srun --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
     --account ibp_pipeline_cleansumstats --pty /bin/bash
   ```
3. **Run a specific test** by copying the relevant section

## Available Tests

| Test | Description | Use Case |
|------|-------------|----------|
| TEST 1 | Normal run (full pipeline) | Verify full workflow works |
| TEST 2 | Reduced sumstat | Quick debugging (~30 variants/chr) |
| TEST 3 | Posteriors only (-2 flag) | Skip scoring step |
| TEST 4 | Step-by-step execution | Test modular v2 features |
| TEST 5 | Interactive development | Debugging inside container |
| TEST 6 | Batch processing | Multiple sumstats from list |

## Test Outputs

Test outputs are stored in subdirectories:
- `out_test_normal/` - Full pipeline runs
- `out_test_reduced/` - Reduced sumstat runs
- `out_test_posteriors_only/` - Posteriors-only runs
- `out_test_stepwise/` - Step-by-step runs
- `out_test_dev/` - Interactive development
- `out_test_batch/` - Batch processing

## Tips

- Use **TEST 2** (reduced sumstat) for quick debugging
- Use **TEST 4** (stepwise) to test the new v2 modular features
- Use **TEST 5** (interactive) when you need to debug inside the container
- Add `-d` flag to keep intermediate files for debugging

## Paths to Update

Update these in `test-commands-v2.sh`:
- `pgsfold` - Path to pgscalculator directory
- `indir` - Path to cleansumstats output
- `id` - Sumstat ID to test with (default: 814)

## See Also

- [README-v2.md](../README-v2.md) - Full v2 documentation
- [test-zone-commands-v2.sh](../test-zone-commands-v2.sh) - Batch template script




