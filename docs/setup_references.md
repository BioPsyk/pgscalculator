# References

## LD references
All LD references will be placed in the folder references/ , but when running the tool, it should be possible to point to anywhere in the system, except through symlinks.

### PRSCS
All references can be downloaded from here: https://github.com/getian107/PRScs

```
# Make PRSCS folder
mkdir -p ld-prscs/kgp
mkdir -p ld-prscs/ukb

# Get the European 1000G project (rename for simpler processing)
wget -O ld-prscs/kgp/kgp.tar.gz https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0
cd ld-prscs/kgp
tar -zxvf kgp.tar.gz
rm kgp.tar.gz
cd ldblk_1kg_eur
for file in ldblk_1kg_chr*.hdf5; do
    # Extract everything after "ldblk_1kg_chr"
    newname="${file#ldblk_1kg_chr}"
    
    # Rename the file
    mv "$file" "$newname"
done
cd ..
mv ldblk_1kg_eur/* .
rm -r ldblk_1kg_eur

cd ../..	

# Get the European UKB project
wget -O ld-prscs/ukb/ukb.tar.gz https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0
cd ld-prscs/ukb
tar -zxvf ukb.tar.gz
rm ukb.tar.gz
cd ldblk_ukbb_eur
for file in ldblk_ukbb_chr*.hdf5; do
    # Extract everything after "ldblk_ukbb_chr"
    newname="${file#ldblk_ukbb_chr}"
    
    # Rename the file
    mv "$file" "$newname"
done
cd ..
mv ldblk_ukbb_eur/* .
rm -r ldblk_ukbb_eur
cd ../..	
```

### sbayesr
All references can be downloaded from here: https://github.com/getian107/PRScs

```
# Make sbayesr folder
mkdir -p ld-sbayesr/ukb

wget https://zenodo.org/records/3350914/files/ukbEURu_hm3_sparse.zip
unzip ukbEURu_hm3_sparse.zip



```






# Sumstats
There is a smaller vcf that can be used for smaller testing, but no snps appear to overlap with the join ldref and target genotypes, so I will put a larger sumstat in the reference that can be used for testing instead.

```
mkdir -p references/sumstats
cp -R /home/jesgaaopen/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library_finngen/R8/version_1.6.0/sumstat_FG1970  references/sumstats/

```

