#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include { wf_prscs } from './modules/subworkflow/wf_prscs.nf'
include { wf_sbayesr } from './modules/subworkflow/wf_sbayesr.nf'

// include processes
include { change_build_sumstats } from './modules/process/pr_format_sumstats.nf'


workflow {
 
  if(!params.calc_posterior){
      Channel.fromPath("${params.input}/calc_posteriors/*", type: 'file').view().set { ch_input }
  }else{
    Channel.fromPath("${params.input}/cleaned_GRCh37.gz", type: 'file').set { ch_input_grch37_map }
    Channel.fromPath("${params.input}/cleaned_GRCh38.gz", type: 'file').set { ch_input_grch38 }
 
    // Change sumstat build to GRCh37 metafile
    change_build_sumstats(ch_input_grch38, ch_input_grch37_map)
    change_build_sumstats.out
    .set { ch_input }
  }
    
  // Choose PGS workflow
  if(params.method=="prscs"){wf_prscs(ch_input)
  }else if(params.method=="sbayesr"){wf_sbayesr(ch_input)
  }else{println("method not avail")}
} 

