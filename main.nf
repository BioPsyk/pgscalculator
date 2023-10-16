#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include { wf_prscs } from './modules/subworkflow/wf_prscs.nf'

// support files from assets
if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }

//"/faststorage/jail/project/proto_psych_pgs/data/reference_gwas/PGC_SCZ_2014.vcf.gz"
//input = file(params.input, checkIfExists: true)
Channel.fromPath(params.input, type: 'file').set { ch_input }

workflow {
    if(params.method=="prscs"){wf_prscs(ch_input, mapfile)}
} 

//params.method
//params.trait
//params.split_gwas_path
//params.n





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



