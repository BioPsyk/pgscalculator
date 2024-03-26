#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 variant_map_for_sbayesr
} from '../process/pr_variant_map_calculations.nf'

workflow wf_map_source_variants {

  take:
  input
  
  main:

  input.view()

}
