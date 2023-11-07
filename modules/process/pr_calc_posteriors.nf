// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'big_mem'

    input:
        tuple val(chr), path(gwas), path(lddir), path("geno.bed"), path("geno.bim"), path("geno.fam")
            val(N)
            
    output:
        tuple val(chr), path("prscs_pst_eff_a1_b0.5_phiauto_chr${chr}.txt")

    script:
        def integerN = Math.round(N as float)
        """
        mkdir prscs

        # Count the number of lines in the file
        num_lines=\$(head -n20 "$gwas" | wc -l)

        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then

          python /repos/PRScs/PRScs.py \
            --ref_dir=$lddir \
            --sst_file=${gwas} \
            --bim_prefix="geno" \
            --n_gwas=${integerN} \
            --chrom=${chr} \
            --out_dir=prscs
        else
          touch prscs_pst_eff_a1_b0.5_phiauto_chr${chr}.txt
        fi
        """ 
}

// calculate per chromosome posterior SNP effects for sBayesR
process calc_posteriors_sbayesr {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev

    label 'big_mem'
    cpus 6
    
    input:
        tuple val(chr), path(gwas_chr), path(ldbin), path(ldinfo)
    
    output:
        tuple val(chr), path("chr${chr}.snpRes")

    script:
        ld_prefix="band_chr${chr}.ldm.sparse"
        """
        # Count the number of lines in the file
        num_lines=\$(head -n20 "$gwas_chr" | wc -l)
        
        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          gctb --sbayes R \
            --gwas-summary ${gwas_chr} \
            --ldm ${ld_prefix} \
            --gamma 0.0,0.01,0.1,1 \
            --pi 0.95,0.02,0.02,0.01 \
            --burn-in 2000 \
            --chain-length 10000 \
            --out-freq 10 \
            --unscale-genotype \
            --exclude-mhc \
            --p-value 0.99 \
            --rsq 0.95 \
            --thread 6 \
            --seed 80851 \
            --out chr${chr}
        else
          touch "chr${chr}.snpRes"
        fi
        """
}


