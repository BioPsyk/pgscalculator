#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_sumstats } from '../process/pr_format_sumstats.nf'

workflow wf_prscs {

    take:
    input
    mapfile
    trait
    n

    main:
    format_sumstats(input, mapfile, "prscs")
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[0], parts[1].replace(".tsv", "")]
    }
    .view()
    .set { result }

    emit:
    result
}


