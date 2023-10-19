# ibp_pgs_pipelines_2



```
# quick start (run prs-cs)
./pgscalculator.sh \
  -s tests/example_data/sumstats/cleaned_GRCh37.gz \
  -l references/ld-prscs/kgp \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/prscs_dir_mapfiles.txt \
  -m "prscs" \
  -c conf/prscs_default_20231019.config \
  -o out

```

