#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_sumstats } from '../process/pr_format_sumstats.nf'
include { calc_posteriors_prscs } from '../process/pr_calc_posteriors.nf'
include { calc_score } from '../process/pr_calc_score.nf'

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

    // channel of genotypes
    Channel.fromPath("${genofile}")
    .splitCsv(sep: '\t', header: false)
    .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
    .map { row -> tuple(row[0], row[1], file("$genodir/${row[2]}")) }
    .filter { type -> type[1] in ['bed', 'bim', 'fam'] }
    .groupTuple()
    .map { chrid, _, files -> [chrid, *files] }
    .set { genotypes }

    // format sumstat
    format_sumstats(input, mapfile, "prscs")
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .set { sumstats }

    //set ld
    Channel.fromPath("${lddir}")
    .set { ch_lddir }

    // Calc posteriors
    sumstats
    .combine(ch_lddir)
    .join(genotypes)
    .set{ ch_calc_posterior_input }
    calc_posteriors_prscs(ch_calc_posterior_input, n)

    
    // Calc score
    calc_posteriors_prscs.out
    .join(genotypes)
    .set{ ch_calc_score_input }
    calc_score(ch_calc_score_input, "${params.prscs_posterior_columns}")

    
    emit:
    calc_score.out
}


