#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { 
 format_sumstats 
 split_on_chromosome
} from '../process/pr_format_sumstats.nf'
include { split_on_chromosome_prscs_ld } from '../process/pr_format_ld.nf'
include { calc_posteriors_prscs } from '../process/pr_calc_posteriors.nf'
include { calc_score } from '../process/pr_calc_score.nf'
include { extract_metadata_from_sumstat } from '../process/pr_extract_metadata_from_sumstat.nf'
include {  
 variant_map_for_prscs
 filter_sumstat_variants_on_map_file
} from '../process/pr_variant_map_calculations.nf'

workflow wf_prscs_calc_posteriors {

    take:
    input

    main:

    // support files from assets
    if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }

    // Make quick access metafile object
//    extract_metadata_from_sumstat(mapfile)

    //set ld dir and ld-metafile
    Channel.fromPath("${params.lddir}")
    .set { ch_lddir }
    Channel.fromPath("${params.lddir}/snpinfo_1kg_hm3")
    .set { ch_lddir_meta }

    // set ld split
    split_on_chromosome_prscs_ld(ch_lddir_meta)
    split_on_chromosome_prscs_ld.out
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .set { ch_lddir_split }

    // channel of genotype bim files
    Channel.fromPath("${params.genofile}")
    .splitCsv(sep: '\t', header: false)
    .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
    .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
    .filter { type -> type[1] in ['bim'] }
    .groupTuple()
    .map { chrid, _, files -> [chrid, *files] }
    .set { genotypes_bim }

    // Split sumstat per chromosome
    split_on_chromosome(input)
    split_on_chromosome.out
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .set { ch_split }

     //ch_split.view()
     //genotypes_bim.view()

    // join for variant map
    ch_split
    .join(genotypes_bim)
    .join(ch_lddir_split)
    .set { ch_to_map }  

    // make variant map
    variant_map_for_prscs(ch_to_map)

    //Filter sumstat based on map
    ch_split
    .join(variant_map_for_prscs.out)
    .set { to_sumstat_variant_filter }
    filter_sumstat_variants_on_map_file(to_sumstat_variant_filter)

 //   // format sumstat
 //   format_sumstats(input, mapfile, "prscs")
 //   .flatMap { it }
 //   .map { file ->
 //     def parts = file.name.split("_")
 //     [parts[1].replace(".tsv", ""), file]
 //   }
 //   .set { sumstats }


 //   // Calc posteriors
 //   sumstats
 //   .combine(ch_lddir)
 //   .join(genotypes)
 //   .combine(extract_metadata_from_sumstat.out)
 //   .set{ ch_calc_posterior_input }
 //   calc_posteriors_prscs(ch_calc_posterior_input)

 //   calc_posteriors_prscs.out.set { ch_calculated_posteriors }
 //   
 //   emit:
 //   ch_calculated_posteriors
    
}

workflow wf_prscs_calc_score {

    take:
    input

    main:

    // channel of genotypes
    Channel.fromPath("${params.genofile}")
    .splitCsv(sep: '\t', header: false)
    .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
    .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
    .filter { type -> type[1] in ['bed', 'bim', 'fam'] }
    .groupTuple()
    .map { chrid, _, files -> [chrid, *files] }
    .set { genotypes }

    // Calc score
    calc_posteriors_prscs.out
    .join(genotypes)
    .set{ ch_calc_score_input }
    calc_score(ch_calc_score_input, "${params.calc_posteriors_prscs.score_columns}")

    emit:
    calc_score.out
}

//workflow wf_prscs {
//
//    take:
//    input
//    ch_input_metafile
//
//    main:
//
//    // support files from assets
//    if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }
//
//    // Make quick access metafile object
//    extract_metadata_from_sumstat(ch_input_metafile)
//
//    // channel of genotypes
//    Channel.fromPath("${params.genofile}")
//    .splitCsv(sep: '\t', header: false)
//    .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
//    .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
//    .filter { type -> type[1] in ['bed', 'bim', 'fam'] }
//    .groupTuple()
//    .map { chrid, _, files -> [chrid, *files] }
//    .set { genotypes }
//
//    // format sumstat
//    format_sumstats(input, mapfile, "prscs")
//    .flatMap { it }
//    .map { file ->
//      def parts = file.name.split("_")
//      [parts[1].replace(".tsv", ""), file]
//    }
//    .set { sumstats }
//
//    //set ld
//    Channel.fromPath("${params.lddir}")
//    .set { ch_lddir }
//
//    // Calc posteriors
//    sumstats
//    .combine(ch_lddir)
//    .join(genotypes)
//    .combine(extract_metadata_from_sumstat.out)
//    .set{ ch_calc_posterior_input }
//    calc_posteriors_prscs(ch_calc_posterior_input)
//    
//    // Calc score
//    calc_posteriors_prscs.out
//    .join(genotypes)
//    .set{ ch_calc_score_input }
//    calc_score(ch_calc_score_input, "${params.calc_posteriors_prscs.score_columns}")
//    
//    emit:
//    calc_score.out
//}
//
//
