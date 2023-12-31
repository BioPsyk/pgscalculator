#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include { wf_prscs } from './modules/subworkflow/wf_prscs.nf'
include { wf_sbayesr } from './modules/subworkflow/wf_sbayesr.nf'

// include processes
include { change_build_sumstats } from './modules/process/pr_format_sumstats.nf'
include { copyConfigFiles } from './modules/process/pr_details.nf'


workflow {
 
  if(!params.calc_posterior){
      Channel.fromPath("${params.input}/calc_posteriors/*", type: 'file').set { ch_input }
  }else{
    Channel.fromPath("${params.input}/cleaned_GRCh37.gz", type: 'file').set { ch_input_grch37_map }
    Channel.fromPath("${params.input}/cleaned_GRCh38.gz", type: 'file').set { ch_input_grch38 }
  }
    
  // Choose PGS workflow
  if(params.method=="prscs"){wf_prscs(ch_input_grch38, ch_input_grch37_map)
  }else if(params.method=="sbayesr"){wf_sbayesr(ch_input_grch38)
  }else{println("method not avail")}


  // Add run details to output
  copyConfigFiles(file("/pgscalculator/nextflow.config"))


} 

