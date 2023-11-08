# ibp_pgs_pipelines_2


## Quick start
Run using an example file for a subset of chromosome 10

```
./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_1 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -m "sbayesr" \
  -c conf/base.config \
  -o out3 \
  -d

./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_1 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink_old \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -m "sbayesr" \
  -c conf/base.config \
  -o out3 \
  -d

./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_1 \
  -l references/ld-prscs/ldblk_1kg_eur \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -m "prscs" \
  -c conf/base.config \
  -o out4 \
  -d

```

## Run a full size sumstat file

```
# quick start (run prs-cs)
./pgscalculator.sh \
  -i references/sumstats/sumstat_FG1970 \
  -l references/ld-prscs/ldblk_1kg_eur \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -m "prscs" \
  -c conf/base.config \
  -o out \
  -d

# quick start (run sbayesR)
./pgscalculator.sh \
  -i references/sumstats/sumstat_FG1970 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -m "sbayesr" \
  -c conf/base.config \
  -o out2 \
  -d

#  -s references/sumstats/sumstat_5668 \
#  -l /faststorage/project/proto_psych_pgs_devel/data/ukb_50k_bigset_2.8M \
#  -l references/ld-sbayesr/ukb/ukbEURu_hm3_shrunk_sparse
#  -l /faststorage/project/proto_psych_pgs_devel/data/ukbEURu_hm3_shrunk_sparse \
```

## Running on HPC
```
# On GDK start interactive node
srun --mem=40g --ntasks 1 --cpus-per-task 6 --time=9:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
srun --mem=10g --ntasks 1 --cpus-per-task 6 --time=1:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash


```

