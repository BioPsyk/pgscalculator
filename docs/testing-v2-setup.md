# Testing pgscalculator v2.0.0 Setup

This guide helps verify that pgscalculator v2.0.0 is properly set up and can access mounted directories.

## Quick Test: Verify Mount Points

### 1. Start Interactive Container Session

```bash
# Using Singularity (default)
srun --mem=10g --ntasks 1 --cpus-per-task 6 --time=0:30:00 \
  --account ibp_pipeline_cleansumstats \
  --pty /bin/bash

# Then start container
cd /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator
singularity shell \
  --contain \
  --cleanenv \
  -B /faststorage:/faststorage \
  -B /pgscalculator:/pgscalculator \
  sif/ibp-pgscalculator-base_version-2.0.0.sif
```

### 2. Run Mount Test Script

Inside the container:

```bash
# Copy test script into container or run directly
bash /pgscalculator/test-mounts.sh

# Or test manually:
pgscalculator --version
ls -la /pgscalculator/
ls -la /pgscalculator/bin/
```

## Testing with Real Data

### Step 1: Test Prep Steps (Run Once)

```bash
# Start interactive session
srun --mem=10g --ntasks 1 --cpus-per-task 6 --time=1:00:00 \
  --account ibp_pipeline_cleansumstats \
  --pty /bin/bash

# Run prep steps
cd /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator

./pgscalculator-v2.sh \
  -i /faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.6.7/sumstat_TEST_ID \
  -l /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -c /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/conf/sbayesr.config \
  -o /faststorage/project/ibp_pipeline_pgscalculator/test-zone/out_test_v2_prep \
  -g /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/plink \
  -f /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  --steps prep
```

### Step 2: Verify Prep Output

```bash
# Check prep output
ls -la /faststorage/project/ibp_pipeline_pgscalculator/test-zone/out_test_v2_prep/prep/
ls -la /faststorage/project/ibp_pipeline_pgscalculator/test-zone/out_test_v2_prep/prep/whitelist/

# Check config was created
cat /faststorage/project/ibp_pipeline_pgscalculator/test-zone/out_test_v2_prep/config.yaml
```

### Step 3: Test Full Pipeline

```bash
# Run full pipeline for one sumstat
./pgscalculator-v2.sh \
  -i /faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.6.7/sumstat_TEST_ID \
  -l /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -c /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/conf/sbayesr.config \
  -o /faststorage/project/ibp_pipeline_pgscalculator/test-zone/out_test_v2 \
  -g /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/plink \
  -f /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt
```

### Step 4: Test Step-by-Step Execution

```bash
# Run only posteriors (assuming prep already done)
./pgscalculator-v2.sh \
  -i /faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.6.7/sumstat_TEST_ID \
  -l /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -c /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/conf/sbayesr.config \
  -o /faststorage/project/ibp_pipeline_pgscalculator/test-zone/out_test_v2 \
  -g /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/plink \
  -f /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  --steps posteriors --skip-prep
```

## Batch Job Testing

Use the provided `test-zone-commands-v2.sh` script:

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/test-zone
cp /faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/test-zone-commands-v2.sh .

# Edit the script to match your setup, then run:
bash test-zone-commands-v2.sh
```

## Troubleshooting

### Issue: "pgscalculator: command not found"

**Solution**: Ensure the Dockerfile includes:
```dockerfile
COPY bin/ /pgscalculator/bin/
ENV PATH="/pgscalculator/bin:${PATH}"
```

### Issue: "Config file not found"

**Solution**: The wrapper creates `config.yaml` in the output directory. Check:
- Output directory is writable
- Container has write access to mounted output directory

### Issue: "LD directory not found"

**Solution**: Verify:
- LD directory path is correct and accessible
- Directory is properly mounted (check with `ls /pgscalculator/`)

### Issue: "Genotype files not found"

**Solution**: Check:
- Genotype manifest file format matches expected format
- All paths in manifest are relative to genodir
- Files exist and are readable

## Verification Checklist

- [ ] Container image built with v2.0.0 CLI
- [ ] `pgscalculator --version` works inside container
- [ ] All mount points accessible (`/pgscalculator/input`, `/pgscalculator/outdir`, etc.)
- [ ] Config file created successfully
- [ ] Prep steps complete successfully
- [ ] Can run individual steps (prep, posteriors, score)
- [ ] Can skip prep when already completed
- [ ] Output files created in correct locations


