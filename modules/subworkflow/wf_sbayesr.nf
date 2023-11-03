#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { format_sumstats } from '../process/pr_format_sumstats.nf'
include { calc_posteriors_sbayesr } from '../process/pr_calc_posteriors.nf'

workflow wf_sbayesr {

    take:
    input
    mapfile
    lddir

    main:
    // format sumstat
    format_sumstats(input, mapfile, "sbayesr")
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .set { sumstats }

    //channel of ldfiles
    Channel
    .fromPath("${lddir}/*.bin")
    .map { file ->
        // Split the file name by underscores and select the third element
        def chrNumber = file.baseName.split("_")[1].replaceAll(/[^0-9]/, '')
        return tuple(chrNumber, file)
    }
    .set { ldfiles1 }


    Channel
    .fromPath("${lddir}/*.info")
    .map { file ->
        // Split the file name by underscores and select the third element
        def chrNumber = file.baseName.split("_")[1].replaceAll(/[^0-9]/, '')
        return tuple(chrNumber, file)
    }
    .set { ldfiles2 }

    ldfiles1
    .join(ldfiles2)
    .set {ch_ldfiles }

    sumstats
    .join(ch_ldfiles)
    .set { ch_calc_posteriors }  

    //ch_calc_posteriors
    calc_posteriors_sbayesr(ch_calc_posteriors)
    
    emit:
    sumstats

}

