#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score {
    publishDir "${params.outdir}/intermediates/scores", mode: 'rellink', overwrite: true
    label 'low_mem'
    
    input:
        tuple val(method), val(chr), path(snp_posteriors), path("geno.pgen"), path("geno.pvar"), path("geno.psam")

    output:
        tuple val(method), path("chr${chr}.score")
    
    script:
        """
        # Count the number of lines in the file
        num_lines=\$(head -n20 ${snp_posteriors} | wc -l)
        
        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          plink2 --pfile geno \
          --out tmp \
          --memory 1000 \
          --threads 1 \
          --score ${snp_posteriors} 4 2 3 header cols=+scoresums ignore-dup-ids

          awk '{gsub(/^#/, ""); print}' tmp.sscore > chr${chr}.score
        else
          touch chr${chr}.score
        fi
        """
}

process calc_merged_score {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    label 'low_mem'
    
    input:
        tuple val(method), path(chrscores)

    output:
        tuple val(method), path("${method}_raw_score_all.gz")
    
    script:
        """
        echo "FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG	SCORE1_SUM	FILE_SUM" > "${method}_raw_score_all"
        awk -vOFS="	" '
          FNR>1{FILE_SUM[\$1]++; ALLELE_CT[\$1]+=\$3; NAMED_ALLELE_DOSAGE_SUM[\$1]+=\$4; SCORE1_SUM[\$1]+=\$6}
          END{for(k in ALLELE_CT){print k, k, ALLELE_CT[k], NAMED_ALLELE_DOSAGE_SUM[k], SCORE1_SUM[k]/ALLELE_CT[k], SCORE1_SUM[k], FILE_SUM[k]}}
        ' ${chrscores} >>  "${method}_raw_score_all"

        gzip "${method}_raw_score_all"

        """
}

