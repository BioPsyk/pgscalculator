// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    label 'big_mem'
    
    input:
        tuple val(chr),
            path(gwas),
            val(N),
            path(ld_bin),
            path(ld_info),
            val(ld_cohort),
            val(ld_population),
            val(ld_snp_set),
            val(ld_format),
            val(bfile),
            path(bed), 
            path(bim), 
            path(fam), 
            val(target_cohort),
            val(target_population),
            val(target_snp_set)
    
    output:
        tuple val(chr),
            path("${target_cohort}_${target_population}_${target_snp_set}_${ld_cohort}_${ld_population}_${ld_snp_set}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt")

    script:
        """
        mkdir ${ld_cohort}_${ld_population}_${ld_snp_set}
        mv ${ld_bin} ${ld_info} ${ld_cohort}_${ld_population}_${ld_snp_set}/
        python ./PRScs.py --ref_dir=\$PWD/${ld_cohort}_${ld_population}_${ld_snp_set} \
            --sst_file=$gwas \
            --bim_prefix=$bfile \
            --n_gwas=$N \
            --chrom=$chr \
            --out_dir=${target_cohort}_${target_population}_${target_snp_set}_${ld_cohort}_${ld_population}_${ld_snp_set}
        """ 
}

