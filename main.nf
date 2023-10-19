#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include { wf_prscs } from './modules/subworkflow/wf_prscs.nf'

// support files from assets
if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }

//"/faststorage/jail/project/proto_psych_pgs/data/reference_gwas/PGC_SCZ_2014.vcf.gz"
//input = file(params.input, checkIfExists: true)
Channel.fromPath(params.input, type: 'file').set { ch_input }

workflow {
  if(params.method=="prscs"){
    wf_prscs(ch_input, 
      mapfile, 
      params.trait, 
      params.n, 
      params.lddir, 
      params.genodir, 
      params.genofile
    )
  }
} 

