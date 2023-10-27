// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    label 'big_mem'
    
    input:
        tuple val(chr), path(gwas), path("geno.bed"), path("geno.bim"), path("geno.fam")
            val(N)
            path(lddir)

//            path(ld_bin),
//            path(ld_info),
//            val(ld_cohort),
//            val(ld_population),
//            val(ld_snp_set),
//            val(ld_format),

//            val(target_cohort),
//            val(target_population),
//            val(target_snp_set)
    
    output:
        tuple val(chr), path("calculated_posterior")

        //path("${target_cohort}_${target_population}_${target_snp_set}_${ld_cohort}_${ld_population}_${ld_snp_set}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt")

    script:
        """

        echo ${lddir}
        mkdir outdir
        python /repos/PRScs/PRScs.py --ref_dir=${lddir} \
            --sst_file=${gwas} \
            --bim_prefix="geno" \
            --n_gwas=${N} \
            --chrom=${chr} \
            --out_dir=outdir
        """ 
}


