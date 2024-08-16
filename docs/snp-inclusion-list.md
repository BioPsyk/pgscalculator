# SNP inclusion list
Can be useful when running on a certain subset of variants, for example a set of variants that have been QCd.

## Using snplist
The file only contains the rsid of interest (picked from the genotype files), for example:
```
rs568182971
rs566854205
rs112920234
rs552582111
rs569167217
rs538292749
rs554788161
rs573863682
rs536478188
rs553300198
```

Apply the snplist Using the -s flag
```
./pgscalculator.sh \
  -j sif/ibp-pgscalculator-base_version-0.5.4.sif \
  -i tests/example_data/sumstats/sumstat_2 \
  -l references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  -g references/genotypes_test/plink \
  -f references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  -c conf/sbayesr.config \
  -s references/genotypes_test/mapfiles/all_snps.txt \
  -o out_test_snpsubset

```

