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
  -s tests/example_data/sumstats/sumstat_1 \
  -l references/ld-prscs/ldblk_1kg_eur \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/prscs_dir_mapfiles.txt \
  -m "sbayesr" \
  -c conf/prscs_default_20231019.config \
  -o out2 \
  -d

```

## Running on HPC:s

```
# On GDK start interactive node
srun --mem=40g --ntasks 1 --cpus-per-task 2 --time=9:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash


```

