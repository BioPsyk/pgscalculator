#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include { wf_prscs } from './modules/subworkflow/wf_prscs.nf'
include { 
  wf_sbayesr_calc_posteriors
  wf_sbayesr_calc_score
 } from './modules/subworkflow/wf_sbayesr.nf'

// include processes
include { change_build_sumstats } from './modules/process/pr_format_sumstats.nf'
include { copyConfigFiles } from './modules/process/pr_details.nf'


workflow {
 
  // Choose PGS workflow
  if(params.method=="prscs"){
    if(params.calc_posterior){
      // Change input sumstat to GRCh37
      Channel.fromPath("${params.input}/cleaned_GRCh37.gz", type: 'file').set { ch_input_grch37_map }
      Channel.fromPath("${params.input}/cleaned_GRCh38.gz", type: 'file').set { ch_input_grch38 }
      change_build_sumstats(ch_input_grch38, ch_input_grch37_map)
      change_build_sumstats.out.set { ch_input }
    }else{
      Channel.fromPath("${params.input}/calc_posteriors/*", type: 'file').set { ch_input }
    }
    wf_prscs(ch_input_grch38, ch_input_grch37_map)
  }else if(params.method=="sbayesr"){
    if(params.calc_posterior){
      Channel.fromPath("${params.input}/cleaned_GRCh38.gz", type: 'file').set { ch_input }
      wf_sbayesr_calc_posteriors(ch_input)
      wf_sbayesr_calc_posteriors.out.ch_calculated_posteriors.set{ ch_calculated_posteriors_2 }
    }else{
      Channel.fromPath("${params.input}/calc_posteriors/*", type: 'file').set { ch_input }
      input
      .map { file ->
        def chrWithPrefix = file.getBaseName().split("_")[0]
        def chr = chrWithPrefix.replaceAll("chr", "")
        return tuple(chr, file)
      }.set { ch_calculated_posteriors_2 }
    }
    if(params.calc_score){
      wf_sbayesr_calc_score(ch_calculated_posteriors_2)
    }
  }else{
    println("method not avail")
  }


  // Add run details to output
  copyConfigFiles(file("/pgscalculator/nextflow.config"))

} 

