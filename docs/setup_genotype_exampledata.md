# 1kgp example genotypes
For development we use 1kgp, which genotypes are not restricted to a secure environmnet.


## Dev software
Use conda for a quick and simple installation of plink

```
mamba create -n prs_dev
mamba install -n prs_dev plink -c bioconda
conda activate prs_dev

```

## Download vcf
```
mkdir -p vcf
cd vcf

for chr in $(seq 1 22); do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz . ;done
for chr in $(seq 1 22); do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi . ;done

cd ..
```

## Convert to plink format bim fam bam

```
mkdir -p plink

for vcf_file in vcf/*.vcf.gz; do
    echo "processing: ${vcf_file}"
    
    # Extract the base name without extension
    base_name=$(basename "$vcf_file" .vcf.gz)
    
    # Convert VCF to PLINK format
    plink --vcf "$vcf_file" --make-bed --out "plink/$base_name"
done

```



