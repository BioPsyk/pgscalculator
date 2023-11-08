#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 add_N_effective
 format_sumstats
 force_EAF_to_sumstat
} from '../process/pr_format_sumstats.nf'
include { calc_posteriors_sbayesr } from '../process/pr_calc_posteriors.nf'
include { calc_score } from '../process/pr_calc_score.nf'

workflow wf_sbayesr {

    take:
    input
    mapfile
    lddir
    metafile
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
    add_N_effective(input, metafile)
    force_EAF_to_sumstat(add_N_effective.out, metafile)
    format_sumstats(force_EAF_to_sumstat.out, mapfile, "sbayesr")
    .flatMap { it }
    .map { file ->
      def parts = file.name.split("_")
      [parts[1].replace(".tsv", ""), file]
    }
    .set { sumstats }

    //channel of ldfiles
    //Channel.fromPath("${lddir}/*.ukb10k.mldm")
    //.set { ch_mldm }

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

    // Calc score
    calc_posteriors_sbayesr.out
    .join(genotypes)
    .set{ ch_calc_score_input }
    calc_score(ch_calc_score_input, "${params.calc_posteriors_sbayesr.score_columns}")
    
    emit:
    calc_score.out

}

