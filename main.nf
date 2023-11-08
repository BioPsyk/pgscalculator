#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include { wf_prscs } from './modules/subworkflow/wf_prscs.nf'
include { wf_sbayesr } from './modules/subworkflow/wf_sbayesr.nf'

// include processes
include { extract_metadata_from_sumstat } from './modules/process/pr_extract_metadata_from_sumstat.nf'
include { change_build_sumstats } from './modules/process/pr_format_sumstats.nf'

// support files from assets
if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }

//input = file(params.input, checkIfExists: true)
Channel.fromPath("${params.input}/cleaned_GRCh37.gz", type: 'file').set { ch_input_grch37_map }
Channel.fromPath("${params.input}/cleaned_GRCh38.gz", type: 'file').set { ch_input_grch38 }
Channel.fromPath("${params.input}/cleaned_metadata.yaml", type: 'file').set { ch_input_metafile }

workflow {
    
  // Change build to GRCh37 metafile
  change_build_sumstats(ch_input_grch38, ch_input_grch37_map)
  
  // Make quick access metafile object
  extract_metadata_from_sumstat(ch_input_metafile)

  // Choose PGS workflow
  if(params.method=="prscs"){
    wf_prscs(change_build_sumstats.out, 
      mapfile, 
      extract_metadata_from_sumstat.out,
      params.lddir, 
      params.genodir, 
      params.genofile
    )
  }else if(params.method=="sbayesr"){
    wf_sbayesr(change_build_sumstats.out, 
      mapfile,
      params.lddir,
      ch_input_metafile,
      params.genodir, 
      params.genofile
    )

  }
} 

