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

wget https://cnsgenomics.com/data/GCTB/band_ukb_10k_hm3.zip
unzip band_ukb_10k_hm3.zip


```






## Sumstats
There is a smaller vcf that can be used for smaller testing, but no snps appear to overlap with the join ldref and target genotypes, so I will put a larger sumstat in the reference that can be used for testing instead.

```
mkdir -p references/sumstats
cp -R /home/jesgaaopen/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library_finngen/R8/version_1.6.0/sumstat_FG1970  references/sumstats/

```

## Genotypes

### rsid mapping
We are using the dbsnp pre-sorted mapfiles from cleansumstats, but compressed using gzip. We have also prepared a sbayesr ld reference specific map, which will be much faster to map against, and only keep the genotypes for which rsids that actually will be used.

```
# Copy from cleansumstats pipeline and place here:
gzip -c <location>/All_20180418_GRCh37_GRCh38.sorted.bed > references/cleansumstat_rsid_map/All_20180418_GRCh37_GRCh38.sorted.bed.gz
gzip -c <location>/All_20180418_GRCh38_GRCh37.sorted.bed > references/cleansumstat_rsid_map/All_20180418_GRCh38_GRCh37.sorted.bed.gz

# Extract list of rsids from the sbayesr reference (FILE 1)
rsidall="references/cleansumstat_rsid_map/All_20180418_GRCh37_GRCh38.sorted.bed.gz"
outfold=references/cleansumstat_rsid_map/sbayesr_band_ukb_10k_hm3/b37
mkdir -p ${outfold}
for chr in $(seq 1 22); do
   (
      echo "Processing chr${chr}"
      date
      ldfile="references/ld-sbayesr/ukb/band_ukb_10k_hm3/band_chr${chr}.ldm.sparse.info"
      rsidchr="${outfold}/chr${chr}_GRCh37_GRCh38.sorted.bed.gz"
      awk '
        NR==FNR && NR>1{a[$2];next}
        NR!=FNR{if($3 in a){print $0}}
      ' "${ldfile}" <(zcat "${rsidall}") | gzip -c > "${rsidchr}"
   ) &
done
wait # This will wait for all background processes to finish

# Extract list of rsids from the sbayesr reference (FILE 1)
rsidall="references/cleansumstat_rsid_map/All_20180418_GRCh38_GRCh37.sorted.bed.gz"
outfold=references/cleansumstat_rsid_map/sbayesr_band_ukb_10k_hm3/b38
mkdir -p ${outfold}
for chr in $(seq 1 22); do
   (
      echo "Processing chr${chr}"
      date
      ldfile="references/ld-sbayesr/ukb/band_ukb_10k_hm3/band_chr${chr}.ldm.sparse.info"
      rsidchr="${outfold}/chr${chr}_GRCh38_GRCh37.sorted.bed.gz"
      awk '
        NR==FNR && NR>1{a[$2];next}
        NR!=FNR{if($3 in a){print $0}}
      ' "${ldfile}" <(zcat "${rsidall}") | gzip -c > "${rsidchr}"
   ) &
done
wait # This will wait for all background processes to finish



```

