# pgscalculator

## Set up

```
# i. Make sure git and singularity are installed, see [singularity installation](docs/singularity-installation.md)
singularity --version
git --version

# ii. clone and enter the cleansumstats github project
git clone https://github.com/BioPsyk/pgscalculator.git
cd pgscalculator

# iii. Download our container image, move it to a folder called tmp within the repo (<1GB)
mkdir -p tmp
chmod ug+rwX tmp
singularity pull tmp/ibp-pgscalculator-base_version-0.5.4.sif docker://biopsyk/ibp-pgscalculator:0.5.4

```

## Quick start

Run using an example file for a random subset of variants on all chromosomes. Just replace the input in -i with a folder cleaned by cleansumstats.

```
# Run both calc posterior and score in one run
./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out3 \
  -d

```

## Running on HPC
examles of running an interactive node on a slurm cluster
```
# On GDK start interactive node (minimum 6 cpus 10g)
srun --mem=10g --ntasks 1 --cpus-per-task 6 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
srun --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash

```

