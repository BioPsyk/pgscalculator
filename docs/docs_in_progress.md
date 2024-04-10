Change config to prscs method
```
# Run both calc posterior and score in one run
./pgscalculator.sh \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-prscs/ldblk_1kg_eur \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/prscs.config \
  -o out3b \
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
  -o out4 \
  -2

# Run only calc score (-1) (-l not required)
# -i is not pointing to the output folder of run only calc posterior
./pgscalculator.sh \
  -i out4 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -o out5 \
  -1

```
