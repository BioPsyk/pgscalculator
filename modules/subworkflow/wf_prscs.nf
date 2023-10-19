#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_sumstats } from '../process/pr_format_sumstats.nf'

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
    // channel of ldfiles
    Channel
    .fromPath("${lddir}/*.hdf5")
    .map { file -> 
        def chrNumber = file.baseName
        return tuple(chrNumber, file) 
    }
    .view()

    // channel of genotypes
    Channel.fromPath("${genofile}")
    .splitCsv(sep: '\t', header: false)
    .map { row -> tuple(row[0], row[1], file(row[2])) }
    .view()

  
    // format sumstat
    format_sumstats(input, mapfile, "prscs")
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .view()
    .set { result }

    emit:
    result
}


