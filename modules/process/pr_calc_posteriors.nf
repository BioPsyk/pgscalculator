// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'big_mem'

    input:
        tuple val(chr), path(gwas), path(lddir), path("geno.bed"), path("geno.bim"), path("geno.fam")
            val(N)
            
    output:
        tuple val(chr), path("outdir_pst_eff_a1_b0.5_phiauto_chr${chr}.txt")

    script:
        def integerN = Math.round(N as float)
        """
        mkdir outdir
        python /repos/PRScs/PRScs.py \
            --ref_dir=$lddir \
            --sst_file=${gwas} \
            --bim_prefix="geno" \
            --n_gwas=${integerN} \
            --chrom=${chr} \
            --out_dir=outdir
        """ 
}


