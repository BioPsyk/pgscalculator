#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score {
    publishDir "${params.outdir}/intermediates/scores", mode: 'rellink', overwrite: true
    label 'mod_mem'
    
    input:
        tuple val(chr), path(snp_posteriors), path("geno.pgen"), path("geno.pvar"), path("geno.psam")

    output:
        path("chr${chr}.score")
    
    script:
        """
        # Count the number of lines in the file
        num_lines=\$(head -n20 ${snp_posteriors} | wc -l)
        
        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          plink2 --pfile geno \
          --out tmp \
          --threads 1 \
          --score ${snp_posteriors} 4 2 3 header cols=+scoresums ignore-dup-ids

          awk '{gsub(/^#/, ""); print}' tmp.sscore > chr${chr}.score
        else
          touch chr${chr}.score
        fi
        """
}

process calc_merged_score {
    publishDir "${params.outdir}/extra", mode: 'copy', overwrite: true
    label 'mod_mem'
    
    input:
        path(chrscores)

    output:
        path("raw_score_all")
    
    script:
        """
        echo "FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG	SCORE1_SUM	FILE_SUM" > raw_score_all
        awk -vOFS="	" '
          FNR>1{FILE_SUM[\$1]++; ALLELE_CT[\$1]+=\$3; NAMED_ALLELE_DOSAGE_SUM[\$1]+=\$4; SCORE1_SUM[\$1]+=\$6}
          END{for(k in ALLELE_CT){print k, k, ALLELE_CT[k], NAMED_ALLELE_DOSAGE_SUM[k], SCORE1_SUM[k]/ALLELE_CT[k], SCORE1_SUM[k], FILE_SUM[k]}}
        ' ${chrscores} >>  "raw_score_all"

        """
}

// rename to indep_pairwise_filter_for_benchmark
process indep_pairwise_for_benchmark {
    publishDir "${params.outdir}/intermediates/indep_pairwise_for_benchmark", mode: 'rellink', overwrite: true
    
    input:
        tuple val(chr), path(sumstat), path("geno.pgen"), path("geno.pvar"), path("geno.psam")

    output:
        tuple val(chr), path("chr${chr}_sumstat")
    
    script:
        """
        awk 'NR>1{print \$4}' ${sumstat} > snps
        plink2 --pfile geno \
         --indep-pairwise 500kb 1 0.2 \
         --extract snps \
         --out "chr${chr}"

        awk 'NR==FNR{a[\$1]; next} FNR==1 || (\$4 in a)' "chr${chr}.prune.in" ${sumstat} > "chr${chr}_sumstat" 
        """
}

