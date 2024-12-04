# Docker
Run using docker image (use the tag: dockerhub_biopsyk)
```bash
## pull docker image
docker pull biopsyk/ibp-pgscalculator:0.5.4

## Run using local docker build 
./pgscalculator.sh \
  -j docker \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out_test_5 

## Run using dockerhub image
./pgscalculator.sh \
  -j dockerhub_biopsyk \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out_test_5 

```
