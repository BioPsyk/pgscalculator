#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include subworkflows
include {
  wf_map_source_variants
 } from './modules/subworkflow/wf_map_source_variants.nf'
include {
  wf_prscs_calc_posteriors
  wf_prscs_calc_score
 } from './modules/subworkflow/wf_prscs.nf'
include { 
  wf_sbayesr_calc_posteriors
  wf_sbayesr_calc_score
 } from './modules/subworkflow/wf_sbayesr.nf'
include { 
  wf_sbayesrc_calc_posteriors
  wf_sbayesrc_calc_score
 } from './modules/subworkflow/wf_sbayesrc.nf'

// include processes
include { add_build_sumstats } from './modules/process/pr_format_sumstats.nf'
include { copyConfigFiles } from './modules/process/pr_details.nf'
include { copyVersionFile } from './modules/process/pr_details.nf'

workflow {

  // Add b37 to sumstat regardless if it is used by the PGS method or not
  Channel.fromPath("${params.input}/cleaned_GRCh38.gz", type: 'file').set { ch_input_grch38 }
  Channel.fromPath("${params.input}/cleaned_GRCh37.gz", type: 'file').set { ch_input_grch37_map }
  add_build_sumstats(ch_input_grch38, ch_input_grch37_map)
  add_build_sumstats.out.set { ch_input }
 
  // Choose PGS workflow
  if(params.method=="prscs"){

    if(!params.calc_posterior){
      Channel.fromPath("${params.input}/calc_posteriors/*", type: 'file').set { ch_input }
      input
      .map { file ->
        def chrWithPrefix = file.getBaseName().split("_")[0]
        def chr = chrWithPrefix.replaceAll("chr", "")
        return tuple(chr, file)
      }.set { ch_calculated_posteriors }
    }else{
      wf_prscs_calc_posteriors(ch_input)
      wf_prscs_calc_posteriors.out.ch_calculated_posteriors.set{ ch_calculated_posteriors }
    }
    if(params.calc_score){
      wf_prscs_calc_score(ch_calculated_posteriors)
    }

  }else if(params.method=="sbayesr"){

    if(!params.calc_posterior){
      Channel.fromPath("${params.input}/calc_posteriors/*", type: 'file').set { ch_input }
      input
      .map { file ->
        def chrWithPrefix = file.getBaseName().split("_")[0]
        def chr = chrWithPrefix.replaceAll("chr", "")
        return tuple(chr, file)
      }.set { ch_calculated_posteriors }
    }else{
      wf_sbayesr_calc_posteriors(ch_input)
      wf_sbayesr_calc_posteriors.out.ch_formatted_posteriors.set{ ch_formatted_posteriors }
      wf_sbayesr_calc_posteriors.out.variant_maps_for_sbayesr.set{ ch_variant_maps }
      wf_sbayesr_calc_posteriors.out.sumstats_filtered.set{ ch_sumstat }
      wf_sbayesr_calc_posteriors.out.unified_plink2_files.set{ ch_unified_plink2 }
    }

    if(params.calc_score){
      wf_sbayesr_calc_score(ch_formatted_posteriors, ch_variant_maps, ch_sumstat, ch_unified_plink2)
    }

  }else if(params.method=="sbayesrc"){

    if(!params.calc_posterior){
      Channel.fromPath("${params.input}/calc_posteriors/*", type: 'file').set { ch_input }
      input
      .map { file ->
        def chrWithPrefix = file.getBaseName().split("_")[0]
        def chr = chrWithPrefix.replaceAll("chr", "")
        return tuple(chr, file)
      }.set { ch_calculated_posteriors }
    }else{
      wf_sbayesrc_calc_posteriors(ch_input)
      wf_sbayesrc_calc_posteriors.out.ch_formatted_posteriors.set{ ch_formatted_posteriors }
      wf_sbayesrc_calc_posteriors.out.variant_maps_for_sbayesr.set{ ch_variant_maps }
      wf_sbayesrc_calc_posteriors.out.sumstats_filtered.set{ ch_sumstat }
      wf_sbayesrc_calc_posteriors.out.unified_plink2_files.set{ ch_unified_plink2 }
    }

    if(params.calc_score){
      wf_sbayesrc_calc_score(ch_formatted_posteriors, ch_variant_maps, ch_sumstat, ch_unified_plink2)
    }

  }else{
    println("method not avail")
  }


  // Add run details to output
  copyConfigFiles(file("/pgscalculator/nextflow.config"))
  copyVersionFile(file("/pgscalculator/VERSION"))

} 

