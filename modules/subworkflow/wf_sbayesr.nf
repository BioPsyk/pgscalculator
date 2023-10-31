#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_sumstats } from '../process/pr_format_sumstats.nf'

workflow wf_sbayesr {

    take:
    input
    mapfile

    main:
    // format sumstat
    format_sumstats(input, mapfile, "sbayesr")
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .set { sumstats }

    // channel of ldfiles
    Channel
    .fromPath("${lddir}/*.hdf5")
    .map { file -> 
        def chrNumber = file.baseName
        return tuple(chrNumber, file) 
    }
    .set { ldfiles }
    

    
    emit:
    sumstats

}


