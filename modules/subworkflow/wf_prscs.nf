#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_sumstats } from '../process/pr_format_sumstats.nf'
include { calc_posteriors_prscs } from '../process/pr_calc_posteriors.nf'

workflow wf_prscs {

    take:
    input
    mapfile
    trait
    n
    lddir
    genodir
    genofile

    main:
   // // channel of ldfiles
   // Channel
   // .fromPath("${lddir}/*.hdf5")
   // .map { file -> 
   //     def chrNumber = file.baseName
   //     return tuple(chrNumber, file) 
   // }
   // .set { ldfiles }

    // channel of genotypes
    Channel.fromPath("${genofile}")
    .splitCsv(sep: '\t', header: false)
    .map { row -> tuple(row[0], row[1], file("$genodir/${row[2]}")) }
    .filter { type -> type[1] in ['bed', 'bim', 'fam'] }
    .groupTuple()
    .map { chrid, _, files -> [chrid, *files] }
    .set { genotypes }

    "${genodir}/snpinfo_1kg_hm3"

    // format sumstat
    format_sumstats(input, mapfile, "prscs")
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .set { sumstats }

    // Calc posteriors
    
    sumstats
    .join(genotypes)
    .set{ calc_posterior_input }
 
    calc_posteriors_prscs(calc_posterior_input, n, file("${lddir}")) \

    //Run PRS-CS

//
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

    

    emit:
    calc_posterior_input
}


