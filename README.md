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
singularity pull sif/ibp-pgscalculator-base_version-0.5.4.sif docker://biopsyk/ibp-pgscalculator:0.5.4
# Run both calc posterior and score in one run
./pgscalculator.sh \
  -j sif/ibp-pgscalculator-base_version-0.5.4.sif \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out_test_1 \
  -d

```


Change config to prscs method
```
# Run both calc posterior and score in one run
./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-prscs/ldblk_1kg_eur \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/prscs.config \
  -o out_test_2 \
  -d

```
## Divide into two 

```
# Run only calc posterior (-2)
./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_2 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -c conf/sbayesr.config \
  -o out_test_3 \
  -2

# Run only calc score (-1) (-l not required)
# -i is not pointing to the output folder of run only calc posterior
./pgscalculator.sh \
  -i out4 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out_test_4 \
  -1

```

### Docker
using docker image (use the tag: dockerhub_biopsyk)
```bash
## pull docker image
docker pull biopsyk/ibp-pgscalculator:0.5.4

## Run using docker tag
./pgscalculator.sh \
  -j dockerhub_biopsyk \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out_test_5 

```
# On GDK start interactive node (minimum 6 cpus 10g)
srun --mem=10g --ntasks 1 --cpus-per-task 6 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
srun --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash

```

