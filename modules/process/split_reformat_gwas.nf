#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Nextflow process to split a gwas cleansumstats default output into chunks
// Generates a space delimited file for prs-cs input
// Generates a tab delimited format format for prs-cs


process split_reformat_gwas {
  label 'mod_mem'

  input:
    tuple path(input_file.tsv),
      val(method),
      path(split_gwas)

  output:
    path('formatted_*')

  script:
    """
    format_sumstats.sh input_file.tsv map_file.csv ${method} formatted_
    """
}

