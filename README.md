# ibp_pgs_pipelines_2


## Quick start
Run using an example file for a subset of chromosome 10. Just replace the input in -i with a folder cleaned by cleansumstats.

```
# Run both calc posterior and score in one run
./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_1 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink_old \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out3

```
## Divide into two 

```
# Run only calc posterior (-2)
./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_1 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -c conf/sbayesr.config \
  -o out4 \
  -2

# Run only calc score (-1)
# -i is not pointing to the output folder of run only calc posterior
./pgscalculator.sh \
  -i out4 \
  -g references/genotypes_test/plink_old \
  -b "37" \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out4 \
  -1

```

## Running on HPC
```
# On GDK start interactive node (minimum 6 cpus 10g)
srun --mem=10g --ntasks 1 --cpus-per-task 6 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash

```

