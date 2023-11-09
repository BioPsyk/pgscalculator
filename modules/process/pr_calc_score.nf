#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score {
    publishDir "${params.outdir}/scores", mode: 'rellink', overwrite: true
    label 'mod_mem'
    
    input:
        tuple val(chr), path(snp_posteriors), path("geno.bed"), path("geno.bim"), path("geno.fam")
        val(snp_posteriors_cols)

    output:
        path("chr${chr}.score")
    
    script:
        """
        # Count the number of lines in the file
        num_lines=\$(head -n20 ${snp_posteriors} | wc -l)
        
        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          plink2 --bfile geno \
          --out tmp \
          --score ${snp_posteriors} ${snp_posteriors_cols} header cols=+scoresums ignore-dup-ids

          awk '{gsub(/^#/, ""); print}' tmp.sscore > chr${chr}.score
        else
          touch chr${chr}.score
        fi
        """
}


