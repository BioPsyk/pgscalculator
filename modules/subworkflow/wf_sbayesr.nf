#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 variant_map_for_sbayesr
 filter_sumstat_variants_on_map_file
} from '../process/pr_variant_map_calculations.nf'

include {  
 rmcol_build_sumstats
 add_N_effective
 format_sumstats
 force_EAF_to_sumstat
 add_B_and_SE
 filter_bad_values_1
 filter_bad_values_2
 filter_on_ldref_rsids
 split_on_chromosome
 concatenate_sumstat_input
} from '../process/pr_format_sumstats.nf'
include { 
 calc_posteriors_sbayesr 
 concatenate_sbayes_posteriors
 qc_posteriors
 format_sbayesr_posteriors
} from '../process/pr_calc_posteriors.nf'
include {
  calc_score 
  calc_merged_score
} from '../process/pr_calc_score.nf'
include { 
  add_rsid_to_genotypes 
  concat_genotypes
} from '../process/pr_format_genotypes.nf'
include { 
  extract_maf_from_genotypes
  concatenate_plink_maf
} from '../process/pr_extract_from_genotypes.nf'

workflow wf_sbayesr_calc_posteriors {

  take:
  input

  main:

  // Read ld from reference
  Channel
  .fromPath("${params.lddir}/*.bin")
  .map { file ->
      // Split the file name by underscores and select the third element
      def chrNumber = file.baseName.split("_")[1].replaceAll(/[^0-9]/, '')
      return tuple(chrNumber, file)
  }
  .set { ldfiles1 }
  Channel
  .fromPath("${params.lddir}/*.info")
  .map { file ->
      // Split the file name by underscores and select the third element
      def chrNumber = file.baseName.split("_")[1].replaceAll(/[^0-9]/, '')
      return tuple(chrNumber, file)
  }
  .set { ldfiles2 }

  ldfiles1
  .join(ldfiles2)
  .set {ch_ldfiles }

  // Metafile for sumstat
  Channel.fromPath("${params.input}/cleaned_metadata.yaml", type: 'file').set { ch_input_metafile }

  // Support files, default from assets/
  if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }

  // channel of genotype bim files
  Channel.fromPath("${params.genofile}")
  .splitCsv(sep: '\t', header: false)
  .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
  .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
  .filter { type -> type[1] in ['bim'] }
  .groupTuple()
  .map { chrid, _, files -> [chrid, *files] }
  .set { genotypes_bim }


  // Split sumstat per chromosome
  split_on_chromosome(input)
  split_on_chromosome.out
  .flatMap { it }
  .map { file ->
    def parts = file.name.split("_")
    [parts[1].replace(".tsv", ""), file]
  }
  .set { ch_split }

  // add metafile
  ch_split
  .combine(ch_input_metafile)
  .set { ch_split2 }

  // format chr-chunked sumstats
  add_N_effective(ch_split2, "${params.whichN}")
  force_EAF_to_sumstat(add_N_effective.out)
  filter_bad_values_1(force_EAF_to_sumstat.out)
  add_B_and_SE(filter_bad_values_1.out)
  filter_bad_values_2(add_B_and_SE.out)

  // join for variant map
  filter_bad_values_2.out
  .join(genotypes_bim)
  .join(ch_ldfiles)
  .set { ch_to_map }  

  // make variant map
  variant_map_for_sbayesr(ch_to_map)

  //Filter sumstat based on map
  filter_bad_values_2.out
  .join(variant_map_for_sbayesr.out.map)
  .set { to_sumstat_variant_filter }
  filter_sumstat_variants_on_map_file(to_sumstat_variant_filter)
  
  // Remove b38 as it is not needed and will continue to be present in the mapfile
  rmcol_build_sumstats(filter_sumstat_variants_on_map_file.out, 2)

  // Formatting according to sbayesr
  format_sumstats(rmcol_build_sumstats.out, mapfile, "sbayesr")
  format_sumstats.out.set { sumstats }



  sumstats
  .join(ch_ldfiles)
  .set { ch_calc_posteriors }  

  //ch_calc_posteriors
  calc_posteriors_sbayesr(ch_calc_posteriors).set { ch_calculated_posteriors }

  //format and id-mapping for scoring
  ch_calculated_posteriors
  .join(variant_map_for_sbayesr.out.map)
  .set { ch_format_posteriors }  
  format_sbayesr_posteriors(ch_format_posteriors).set { ch_formatted_posteriors }

  // concat all posteriors
  ch_calculated_posteriors.map {x,y -> y}.collect().set {ch_collected_posteriors}
  concatenate_sbayes_posteriors(ch_collected_posteriors)
  concatenate_sbayes_posteriors.out.set { ch_concatenated_posteriors }

  // concat all input (for QC plots)
  filter_bad_values_2.out.map {x,y -> y}.collect().set {ch_collected_input}
  concatenate_sumstat_input(ch_collected_input)
  concatenate_sumstat_input.out.set { ch_concatenated_input }

  // post QC plots
  filter_bad_values_2.out
  .mix(ch_concatenated_input)
  .set { ch_input_for_qc }

  ch_calculated_posteriors
  .mix(ch_concatenated_posteriors)
  .join(ch_input_for_qc, by: 0)
  .set { ch_posteriors_for_qc }

  qc_posteriors(ch_posteriors_for_qc)

  emit:
  ch_formatted_posteriors
  variant_maps_for_sbayesr = variant_map_for_sbayesr.out.map
}

workflow wf_sbayesr_calc_score {

  take:
  ch_formatted_posteriors
  variant_maps_for_sbayesr

  main:

  // channel of genotypes
  Channel.fromPath("${params.genofile}")
  .splitCsv(sep: '\t', header: false)
  .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
  .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
  .filter { type -> type[1] in ['bed', 'bim', 'fam'] }
  .groupTuple()
  .map { chrid, _, files -> [chrid, *files] }
  .set { genotypes }

  genotypes
  .join(variant_maps_for_sbayesr)
  .set { ch_to_extract_maf }

  extract_maf_from_genotypes(ch_to_extract_maf)
  extract_maf_from_genotypes.out.map {x,y -> y}.collect().set {ch_collected_maf}
  concatenate_plink_maf(ch_collected_maf)

  if(params.concat_genotypes){

    // concat all genotypes
    genotypes
    .map { chr, bed, bim, fam -> 
        def canonicalBed = new File(bed.toString()).canonicalPath
        def canonicalBim = new File(bim.toString()).canonicalPath
        def canonicalFam = new File(fam.toString()).canonicalPath
        return "$canonicalBed $canonicalBim $canonicalFam"
    }
    .collectFile(){ content -> [ "allgenotypes.txt", content + '\n' ] }
    .set { ch_all_genofile }
    concat_genotypes(ch_all_genofile)

    // join posteriors and genotypes all chromosomes
    ch_concatenated_posteriors
    .join(concat_genotypes.out)
    .set{ ch_calc_score_input_all }
  }else{
    Channel.empty().set { ch_calc_score_input_all }
  }

  // join posteriors and genotypes per chromosome
  ch_formatted_posteriors
  .join(genotypes)
  .set{ ch_calc_score_input_per_chr }

  // Calc score
  ch_calc_score_input_per_chr
  .mix(ch_calc_score_input_all)
  .set{ ch_calc_score }
  calc_score(ch_calc_score)

  // Merge and calculate per chromosome scores
  calc_merged_score(calc_score.out.collect())
}


