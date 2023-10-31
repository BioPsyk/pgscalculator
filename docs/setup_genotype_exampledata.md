# 1kgp example genotypes
For development we use 1kgp, which genotypes are not restricted to a secure environmnet.


## Dev software
Use conda for a quick and simple installation of plink

```
# reusing this environment for plink and bcftools
mamba create -n merging_imputed_genotypes --channel bioconda \
  nextflow==23.10.0 \
  plink=1.90b6.21 \
  plink2=2.00a5 \
  bcftools=1.18 \
  tabix=1.11

conda activate merging_imputed_genotypes

```

## Download vcf
```
mkdir -p vcf
cd vcf

for chr in $(seq 1 22); do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz . ;done
for chr in $(seq 1 22); do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi . ;done

cd ..
```

## Download 1000g vcf calls with rsids
```
mkdir -p vcf_calls

cd vcf_calls
wget https://ftp.ensembl.org/pub/grch37/release-110/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz
wget https://ftp.ensembl.org/pub/grch37/release-110/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi
cd ..
```

## Annotate vcf with rsids from the vcf calls file
```
mkdir -p vcf_rsid
cd vcf_rsid
conda activate duplicate_id_resolver
sbatch -a 1-22%22 --mem=8g --ntasks 1 --cpus-per-task 2 --time=9:00:00 --account ibp_pipeline_cleansumstats --job-name="mega_%A_%a" --output="mega_%A_%a.out" --error="mega_%A_%a.err" --wrap="
  echo \${SLURM_ARRAY_TASK_ID} 
  date
    bcftools annotate -a ../vcf_calls/1000GENOMES-phase_3.vcf.gz -c ID -Oz ../vcf/ALL.chr\${SLURM_ARRAY_TASK_ID}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz  > ALL.chr\${SLURM_ARRAY_TASK_ID}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
  date
"

cd ..
```

## Convert to plink format bim fam bam
srun --mem=4g --ntasks 1 --cpus-per-task 1 --time=4:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
mkdir -p plink

for vcf_file in vcf_rsid/*.vcf.gz; do
    echo "processing: ${vcf_file}"
    
    # Extract the base name without extension
    base_name=$(basename "$vcf_file" .vcf.gz)
    
    # Convert VCF to PLINK format
    plink --vcf "$vcf_file" --make-bed --out "plink/$base_name"
done

```



