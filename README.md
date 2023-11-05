# ibp_pgs_pipelines_2


## Quick start

```
# quick start (run prs-cs)
./pgscalculator.sh \
  -s references/sumstats/sumstat_FG1970 \
  -l references/ld-prscs/ldblk_1kg_eur \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/prscs_dir_mapfiles.txt \
  -m "prscs" \
  -c conf/prscs_default_20231019.config \
  -o out \
  -d

# quick start (run sbayesR)
./pgscalculator.sh \
  -s references/sumstats/sumstat_5668 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/prscs_dir_mapfiles.txt \
  -m "sbayesr" \
  -c conf/prscs_default_20231019.config \
  -o out2 \
  -d

#  -l /faststorage/project/proto_psych_pgs_devel/data/ukb_50k_bigset_2.8M \
#  -l references/ld-sbayesr/ukb/ukbEURu_hm3_shrunk_sparse
#  -l /faststorage/project/proto_psych_pgs_devel/data/ukbEURu_hm3_shrunk_sparse \
```

## Running on HPC
```
# On GDK start interactive node
srun --mem=40g --ntasks 1 --cpus-per-task 2 --time=9:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash


```

