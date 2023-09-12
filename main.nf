#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { split_reformat_gwas as split_for_prscs } from './modules/split_reformat_gwas.nf'

def help_msg() {
    log.info """
    A nextflow pipeline for performing polygenic score analysis at IBP
    Author: Jesper R. Gådin | jesper.gaadin@gmail.com
    Author: Vivek Appadurai | vivek.appadurai@regionh.dk

    This pipeline is an attempt to simplify the external execution, for example wrapping this nextflow pipeline in a bash wrapper to easier distribute it on different HPCs. Another difference to the previous pipeline is the removal of the methods comparison module, which instead will become its own pipeline using the output from this pipeline.

    Usage: nextflow run main.nf 
    
    Options

    --ref <gwas.vcf.gz> [A reference file of gwas summary stats in VCF.gz format] (Default: PGC SCZ 2014)
    --n <77096> [Reference GWAS Sample Size] (Default: 77096)
    --n_cases <Number of cases> [Number of cases in reference GWAS] (Default: 33640)
    --prevalence <Prevalence in decimals> [If the outcome is binary, population prevalence is used to calculate liability transformed r2] (Default: 0.1)
    --target <iPSYCH2012_All_Imputed_2021_QCed.json> [JSON file of target genotypes to score] (Default: iPSYCH2012 Imputed in 2021)
    --trait <ipsych_scz_2014> [A prefix for output files] (Default: Simple name of reference file)
    --covs <file.covs> [Path to covariates you might want to include in the NULL PGS model, such as age, gender, 10 PCs]
    --pheno <trait.pheno> [Path to the true phenotype file to evaluate the PGS performance]
    --binary <T/F> [Is the outcome binary or continuous?] (Default: T)
    --p_vals <comma separated list of p-values> [thresholds to be use for pruning and thresholding] (Default: 5e-8,1e-6,0.05,1)
    --help <prints this message>
    """
}

params.trait    = file(params.ref).getSimpleName()
target_prefix   = file(params.target).getSimpleName() 

if(params.help)
{
    help_msg()
    exit 0
}

log.info """
============================================================================================================
I B P - P R S -  P I P E L I N E _ v. 1.0 - N F
============================================================================================================
Reference GWAS               : $params.ref
Trait Name                   : $params.trait
Reference GWAS Sample Size   : $params.n
Number of cases              : $params.n_cases
Binary Trait?                : $params.binary
Trait prevalence             : $params.prevalence
Target Genotypes             : $params.target
Target Prefix to use         : $target_prefix
PRSCS 1000G hapmap3 SNPs LD  : $params.prscs_1000G_hm3_eur_ld
PRSCS UKBB hapmap3 LD        : $params.prscs_ukbb_hm3_eur_ld
sBayesR UKBB hapmap3 SNPs LD : $params.sbayesr_ukbb_hm3_eur_ld
sBayesR UKBB 2.5M SNPs LD    : $params.sbayesr_ukbb_big_eur_ld
Output Directory             : $launchDir
Covariates                   : $params.covs
Phenotype                    : $params.pheno
p-value thresholds for P&T   : $params.p_vals
PLINK Path                   : $params.plink_path
Split GWAS Path              : $params.split_gwas_path
sBayesR Path                 : $params.sbayesr_path
PRS_CS Path                  : $params.prscs_path
PRSICE Path                  : $params.prsice_path
PRSICE column checker        : $params.prsice_col_checker
============================================================================================================
"""

ref_ch = Channel.of(1..22) 
    | map {a -> [a, params.ref, "${params.ref}.tbi"]}


workflow {
    Channel.of(1..22) \
    | combine(Channel.of(params.trait)) \
    | combine(ref_ch, by: 0) \
    | combine(Channel.of('prscs')) \
    | combine(Channel.of(params.split_gwas_path)) \
    | combine(Channel.of(params.n)) \
    | split_for_prscs \
    | set { prscs_input_ch } 

    //Run PRS-CS
//
//    Channel.of(1..22) \
//    | combine(prscs_input_ch, by: 0) \
//    | combine(Channel.of(params.n)) \
//    | combine(prscs_ukbb_hm3_eur_ld_ch, by: 0) \
//    | combine(genotypes_ch, by: 0) \
//    | combine(Channel.of(params.trait)) \
//    | combine(Channel.of(params.prscs_path)) \
//    | calc_posteriors_prscs_ukbb_eur_hm3 \
//    | combine(Channel.of("2 4 6")) \
//    | combine(Channel.of("${params.trait}_prscs_ukbb_eur_hm3")) \
//    | combine(genotypes_ch, by: 0) \
//    | combine(Channel.of(params.plink_path)) \
//    | calc_score_prscs_ukbb_eur_hm3 \
//    | collectFile(name: "${target_prefix}_${params.trait}_prscs_ukbb_eur_hm3.score",
//        keepHeader: true,
//        skip: 1) \
//    | set { prscs_ukbb_eur_hm3_score_ch } 
//
//    Channel.of(1..22) \
//    | combine(prscs_input_ch, by: 0) \
//    | combine(Channel.of(params.n)) \
//    | combine(prscs_1kg_hm3_eur_ld_ch, by: 0) \
//    | combine(genotypes_ch, by: 0) \
//    | combine(Channel.of(params.trait)) \
//    | combine(Channel.of(params.prscs_path)) \
//    | calc_posteriors_prscs_1kg_eur_hm3 \
//    | combine(Channel.of("2 4 6")) \
//    | combine(Channel.of("${params.trait}_prscs_1kg_hm3_eur")) \
//    | combine(genotypes_ch, by: 0) \
//    | combine(Channel.of(params.plink_path)) \
//    | calc_score_prscs_1kg_eur_hm3 \
//    | collectFile(name: "${target_prefix}_${params.trait}_prscs_1kg_eur_hm3.score",
//        keepHeader: true,
//        skip: 1) \
//    | set { prscs_1kg_eur_hm3_score_ch }
//

} 


