# pgscalculator

This workflow aims to run different PGS methods using cleansumstats output so the result is comparable.

_Created by Jesper R. Gådin and Andrew Schork_

## Quick start
Make sure git and singularity (or docker) are installed. Then clone the code from github.
```
# Are the required software available 
singularity --version
docker --version
git --version

# Clone and enter the pgscalculator github project
git clone https://github.com/BioPsyk/pgscalculator.git
cd pgscalculator
```

### Singularity
Using singularity (use path to image). Run using an example file for a random subset of variants on all chromosomes. Just replace the input in -i with a folder cleaned by cleansumstats.

```bash
## pull singularity image returning the image as a file (<1GB)
mkdir -p sif
singularity pull sif/ibp-pgscalculator-base_version-0.6.0.sif docker://biopsyk/ibp-pgscalculator:0.6.0-amd64

> **Note**: For ARM64 systems (e.g., Apple Silicon Macs), append `-arm64` to the version tag instead of `-amd64`. For example: `1.3.0-arm64`.

# Run both calc posterior and score in one run
./pgscalculator.sh \
  -j sif/ibp-pgscalculator-base_version-0.5.4.sif \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out_test_1

```

Change config to prscs method
```
# Run both calc posterior and score in one run
./pgscalculator.sh \
  -j sif/ibp-pgscalculator-base_version-0.5.4.sif \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-prscs/ldblk_1kg_eur \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/prscs.config \
  -o out_test_2

```

Change config to sbayesRC method (requires functional annotations)
```bash
# Run both calc posterior and score in one run with sbayesRC
./pgscalculator.sh \
  -j sif/ibp-pgscalculator-base_version-0.6.0.sif \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/sbayesrc/ld-hapmap3/ukbEUR_HM3 \
  -g references/genotypes_test/plink2 \
  -f references/genotypes_test/mapfiles/plink2_genodir_genofiles.txt \
  -c conf/sbayesrc.config \
  -o out_test_3

```

### Specifying Resources
It is possible to run both interactive and batch jobs. Below is an example on GDK(HPC) starting an interactive node. 

**Resource Requirements by Method:**
- **PRS-CS/SBayesR**: Minimum 6 CPUs, 10GB RAM; Recommended 22 CPUs, 20GB RAM
- **SBayesRC**: Minimum 6 CPUs, 20GB RAM; Recommended 22 CPUs, 40GB RAM (higher due to functional annotations)

```bash
# For PRS-CS and SBayesR
srun --mem=10g --ntasks 1 --cpus-per-task 6 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
srun --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash

# For SBayesRC (requires more memory for annotations)
srun --mem=40g --ntasks 1 --cpus-per-task 22 --time=2:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
```

**Note**: SBayesRC runtime varies significantly based on annotation data:
- HapMap3 LD (~1M SNPs): 30 minutes - 2 hours
- Imputed LD (~7M SNPs): 2-8 hours
See [SBayesRC Performance](docs/sbayesrc-performance.md) for detailed benchmarks.

## More documentation
- See [SNP inclusion list](docs/snp-inclusion-list.md)
- See [Two step analysis](docs/two-step-analysis.md)
- See [Use Docker](docs/using-docker.md)
- See [SBayesRC Performance](docs/sbayesrc-performance.md)
- See [FAQ](docs/FAQ.md)


