#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score {
    publishDir "${params.outdir}/scores", mode: 'copy', overwrite: true
    label 'mod_mem'
    
    input:
        tuple val(chr), path(snp_posteriors), path("geno.pgen"), path("geno.pvar"), path("geno.psam")
        val(snp_posteriors_cols)

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
          --score ${snp_posteriors} ${snp_posteriors_cols} header cols=+scoresums ignore-dup-ids

          awk '{gsub(/^#/, ""); print}' tmp.sscore > chr${chr}.score
        else
          touch chr${chr}.score
        fi
        """
}

process calc_merged_score {
    publishDir "${params.outdir}/scores", mode: 'copy', overwrite: true
    label 'mod_mem'
    
    input:
        path(chrscores)

    output:
        path("chrmerged.score")
    
    script:
        """
        echo "FID	IID	ALLELE_CT	NAMED_ALLELE_DOSAGE_SUM	SCORE1_AVG	SCORE1_SUM	FILE_SUM" > chrmerged.score
        awk -vOFS="	" '
          FNR>1{FILE_SUM[\$1]++; ALLELE_CT[\$1]+=\$3; NAMED_ALLELE_DOSAGE_SUM[\$1]+=\$4; SCORE1_SUM[\$1]+=\$6}
          END{for(k in ALLELE_CT){print k, k, ALLELE_CT[k], NAMED_ALLELE_DOSAGE_SUM[k], SCORE1_SUM[k]/ALLELE_CT[k], SCORE1_SUM[k], FILE_SUM[k]}}
        ' ${chrscores} >>  "chrmerged.score"

        """
}

