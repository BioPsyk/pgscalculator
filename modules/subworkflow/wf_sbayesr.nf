#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 variant_map_for_sbayesr
 concatenate_variant_map
 filter_sumstat_variants_on_map_file
 make_snplist_from_bim
 sort_user_snplist
} from '../process/pr_variant_map_calculations.nf'

include {
 rmcol_build_sumstats
 filter_NA_coordinates
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
  make_geno_bim_snpid_unique_bim
  make_geno_bim_snpid_unique_bim_fam_bed
  add_rsid_to_genotypes 
  concat_genotypes
} from '../process/pr_format_genotypes.nf'
include { 
  extract_maf_from_genotypes
  concatenate_plink_maf
} from '../process/pr_extract_from_genotypes.nf'
include { 
  prepare_sumstat_for_benchmark_scoring
  indep_pairwise_for_benchmark
  sumstat_maf_filter
} from '../process/pr_prepare_benchmark.nf'
include { 
  make_augmented_gwas
  concatenate_augmented_sumstat
} from '../process/pr_organize_output.nf'


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
  .set { genotypes_bim_0 }

  make_geno_bim_snpid_unique_bim(genotypes_bim_0)
  make_geno_bim_snpid_unique_bim.out.set { genotypes_bim }

  // SNPlist
  if (params.snplist) {
    Channel.fromPath("${params.snplist}", type: 'file').set { ch_input_snplist_0 }
    sort_user_snplist(ch_input_snplist_0)
    sort_user_snplist.out.set { ch_input_snplist }
  } else {
    make_snplist_from_bim(genotypes_bim.map {x,y -> y}.collect())
    make_snplist_from_bim.out.set { ch_input_snplist }
  }


  // Split sumstat per chromosome
  split_on_chromosome(input)
  split_on_chromosome.out
  .flatMap { it }
  .map { file ->
    def parts = file.name.split("_")
    [parts[1].replace(".tsv", ""), file]
  }
  .filter { item ->
    int chromNumber = item[0].toInteger()
    return chromNumber >= 1 && chromNumber <= 22
  }
  .set { ch_split }


  // Early filter on NA coordinates for b37
  filter_NA_coordinates(ch_split)

  // add metafile
  filter_NA_coordinates.out
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
  .combine(ch_input_snplist)
  .join(ch_ldfiles)
  .set { ch_to_map }  

  // make variant map
  variant_map_for_sbayesr(ch_to_map)
  concatenate_variant_map(variant_map_for_sbayesr.out.map {x,y,z -> y}.collect())

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
  sumstats_filtered=rmcol_build_sumstats.out
  ch_formatted_posteriors
  variant_maps_for_sbayesr = variant_map_for_sbayesr.out.map
}

workflow wf_sbayesr_calc_score {

  take:
  ch_formatted_posteriors
  variant_maps_for_sbayesr
  sumstat

  main:

  // channel of genotypes
  Channel.fromPath("${params.genofile}")
  .splitCsv(sep: '\t', header: false)
  .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
  .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
  .filter { type -> type[1] in ['bed', 'bim', 'fam'] }
  .groupTuple()
  .map { chrid, _, files -> [chrid, *files] }
  .set { genotypes_0 }

  make_geno_bim_snpid_unique_bim_fam_bed(genotypes_0)
  make_geno_bim_snpid_unique_bim_fam_bed.out.set { genotypes }

  // Extract maf from genotypes
  genotypes
  .join(variant_maps_for_sbayesr)
  .set { ch_to_extract_maf }
  extract_maf_from_genotypes(ch_to_extract_maf)
  extract_maf_from_genotypes.out.map {x,y -> y}.collect().set {ch_collected_maf}
  concatenate_plink_maf(ch_collected_maf)

  // prepare benchmark scoring
  sumstat
  .join(variant_maps_for_sbayesr)
  .set { ch_prepare_score_benchmark }
  prepare_sumstat_for_benchmark_scoring(ch_prepare_score_benchmark)

  prepare_sumstat_for_benchmark_scoring.out
  .join(extract_maf_from_genotypes.out)
  .set { ch_prepare_score_benchmark_maf }
  sumstat_maf_filter(ch_prepare_score_benchmark_maf, 0.05)

  sumstat_maf_filter.out
  .join(genotypes)
  .set { ch_for_indep_pairwise }
  indep_pairwise_for_benchmark(ch_for_indep_pairwise)
  indep_pairwise_for_benchmark.out.set { ch_benchmark_ready_to_score }
  
  // not certain this geno concatination will be needed as an option
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

  // Mix main method and benchmark method after adding genotypes
  ch_formatted_posteriors_main = ch_formatted_posteriors.join(genotypes).map { tuple -> ['main'] + tuple }
  ch_benchmark_ready_to_score_bench = ch_benchmark_ready_to_score.join(genotypes).map { tuple -> ['bench'] + tuple }
  ch_calc_score_input_per_chr = ch_formatted_posteriors_main.mix(ch_benchmark_ready_to_score_bench)

  // Calc score
  ch_calc_score_input_per_chr
  .mix(ch_calc_score_input_all)
  .set{ ch_calc_score }
  calc_score(ch_calc_score)

  // Step 1: Branch off and prepare for collection
  ch_main = calc_score.out.filter { it[0] == 'main' }.map { it - 'main' }
  ch_bench = calc_score.out.filter { it[0] == 'bench' }.map { it - 'bench' }
  
  // Step 2: Collect items
  ch_main_collected = ch_main.collect()
  ch_bench_collected = ch_bench.collect()
  
  // Step 3: Process collected items, reintroducing the method info
  ch_main_collected.map { ['main', it] }.set { input_for_main_process }
  ch_bench_collected.map { ['bench', it] }.set { input_for_bench_process }

  // Again Mix main method and benchmark method 
  ch_to_merge_collected = input_for_main_process.mix(input_for_bench_process)

  // Merged chromosomes scores
  calc_merged_score(ch_to_merge_collected)

  // Make Augmented GWAS
  sumstat
  .join(extract_maf_from_genotypes.out)
  .join(ch_formatted_posteriors)
  .join(ch_benchmark_ready_to_score)
  .join(variant_maps_for_sbayesr)
  .set { for_augmented_gwas }
  make_augmented_gwas(for_augmented_gwas)
  concatenate_augmented_sumstat(make_augmented_gwas.out.map {x,y -> y}.collect())

}


