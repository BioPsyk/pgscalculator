#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_sumstats } from '../process/pr_format_sumstats.nf'

workflow wf_prscs {

    take:
    input
    mapfile

    main:
    format_sumstats(input, mapfile, "prscs")
    .set { result }

    emit:
    result
}

