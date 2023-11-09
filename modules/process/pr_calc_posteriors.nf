// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    publishDir "${params.outdir}/calc_posteriors", mode: 'copy', overwrite: true
    label 'big_mem'

    input:
        tuple val(chr), path(gwas), path(lddir), path("geno.bed"), path("geno.bim"), path("geno.fam"), val(N)
            
    output:
        tuple val(chr), path("chr${chr}.posteriors")

    script:
        def integerN = Math.round(N as float)
        """
        mkdir prscs

        # Count the number of lines in the file
        num_lines=\$(head -n20 "$gwas" | wc -l)

        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          # PRScs requires python2, which is best accessed through conda
          #source /opt/micromamba/bin/activate py27

          micromamba run -n py27 python /repos/PRScs/PRScs.py \
            --ref_dir=$lddir \
            --sst_file=${gwas} \
            --bim_prefix="geno" \
            --n_gwas=${integerN} \
            --chrom=${chr} \
            --out_dir=prscs

          mv prscs_pst_eff_a1_b0.5_phiauto_chr${chr}.txt chr${chr}.posteriors
        else
          touch chr${chr}.posteriors
        fi
        """ 
}

// calculate per chromosome posterior SNP effects for sBayesR
process calc_posteriors_sbayesr {
    publishDir "${params.outdir}/calc_posteriors", mode: 'copy', overwrite: true

    label 'big_mem'
    cpus 6
    
    input:
        tuple val(chr), path(gwas_chr), path(ldbin), path(ldinfo)
    
    output:
        tuple val(chr), path("chr${chr}.posteriors")

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
            --gamma ${params.calc_posteriors_sbayesr.gamma} \
            --pi ${params.calc_posteriors_sbayesr.pi} \
            --burn-in ${params.calc_posteriors_sbayesr.burn_in} \
            --chain-length ${params.calc_posteriors_sbayesr.chain_length} \
            --out-freq ${params.calc_posteriors_sbayesr.out_freq} \
            --p-value ${params.calc_posteriors_sbayesr.p_value} \
            --rsq ${params.calc_posteriors_sbayesr.rsq} \
            --thread ${params.calc_posteriors_sbayesr.thread} \
            --seed ${params.calc_posteriors_sbayesr.seed} \
            ${params.calc_posteriors_sbayesr.unscale_genotype} \
            ${params.calc_posteriors_sbayesr.exclude_mhc} \
            --out chr${chr}
          mv "chr${chr}.snpRes" "chr${chr}.posteriors"
        else
          touch "chr${chr}.posteriors"
        fi
        """
}


