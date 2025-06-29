#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 variant_map_for_sbayesr
 concatenate_variant_map
 filter_sumstat_variants_on_map_file
 make_snplist_from_pvar
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
  convert_plink1_to_plink2
  make_geno_pvar_snpid_unique_pvar
  make_geno_pvar_snpid_unique_pvar_psam_pgen
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
  .view { params.dev ? "DEBUG: ld_bin_files_found: $it" : null }
  .map { file ->
      // Split the file name by underscores and select the third element
      def chrNumber = file.baseName.split("_")[1].replaceAll(/[^0-9]/, '')
      return tuple(chrNumber, file)
  }
  .view { params.dev ? "DEBUG: ld_bin_mapped: $it" : null }
  .set { ldfiles1 }
  Channel
  .fromPath("${params.lddir}/*.info")
  .view { params.dev ? "DEBUG: ld_info_files_found: $it" : null }
  .map { file ->
      // Split the file name by underscores and select the third element
      def chrNumber = file.baseName.split("_")[1].replaceAll(/[^0-9]/, '')
      return tuple(chrNumber, file)
  }
  .view { params.dev ? "DEBUG: ld_info_mapped: $it" : null }
  .set { ldfiles2 }

  ldfiles1
  .join(ldfiles2)
  .view { params.dev ? "DEBUG: ch_ldfiles_joined: $it" : null }
  .map { 
    if (params.dev) {
      file("${params.outdir}/channel-trace").mkdirs()
      def chr = it[0]
      def binfile = it[1]
      def infofile = it[2]
      def traceFile = file("${params.outdir}/channel-trace/ch_ldfiles_${chr}.txt")
      traceFile.text = "CHR: ${chr}\nBIN: ${binfile}\nINFO: ${infofile}\n"
    }
    return it
  }
  .set {ch_ldfiles }

  // Metafile for sumstat
  Channel.fromPath("${params.input}/cleaned_metadata.yaml", type: 'file').set { ch_input_metafile }

  // Support files, default from assets/
  if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }

  // Read all genotype files and detect format
  Channel.fromPath("${params.genofile}")
  .splitCsv(sep: '\t', header: false)
  .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
  .view { params.dev ? "DEBUG: genofile_csv_row: $it" : null }
  .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
  .view { params.dev ? "DEBUG: genofile_mapped: $it" : null }
  .groupTuple()
  .view { params.dev ? "DEBUG: genofile_grouped: $it" : null }
  .map { chrid, types, files ->
    // Detect format for this chromosome
    def hasPlink2 = types.any { it in ['pgen', 'pvar', 'psam'] }
    def hasPlink1 = types.any { it in ['bed', 'bim', 'fam'] }
    def format = hasPlink2 ? 'plink2' : (hasPlink1 ? 'plink1' : 'unknown')
    
    if (params.dev) { 
      println("DEBUG: CHR ${chrid} - Format: ${format}, Types: ${types}, Has PLINK2: ${hasPlink2}, Has PLINK1: ${hasPlink1}")
    }
    
    return tuple(chrid, format, types, files)
  }
  .view { params.dev ? "DEBUG: format_detected: $it" : null }
  .branch {
    plink1: it[1] == 'plink1'
    plink2: it[1] == 'plink2'
    unknown: true
  }
  .set { genotype_formats }

  // Handle PLINK1 format - convert to PLINK2
  genotype_formats.plink1
  .map { chrid, format, types, files ->
    // Extract bed, bim, fam files
    def fileMap = [types, files].transpose().collectEntries()
    return tuple(chrid, fileMap['bed'], fileMap['bim'], fileMap['fam'])
  }
  .view { params.dev ? "DEBUG: plink1_to_convert: $it" : null }
  .set { plink1_files }

  convert_plink1_to_plink2(plink1_files)
  
  // Handle native PLINK2 format
  genotype_formats.plink2
  .map { chrid, format, types, files ->
    // Extract pgen, pvar, psam files 
    def fileMap = [types, files].transpose().collectEntries()
    return tuple(chrid, fileMap['pgen'], fileMap['pvar'], fileMap['psam'])
  }
  .view { params.dev ? "DEBUG: native_plink2: $it" : null }
  .set { native_plink2_files }

  // Combine converted and native PLINK2 files
  convert_plink1_to_plink2.out
  .mix(native_plink2_files)
  .view { params.dev ? "DEBUG: unified_plink2_files: $it" : null }
  .map { 
    if (params.dev) {
      file("${params.outdir}/channel-trace").mkdirs()
      def chr = it[0]
      def pgen = it[1]
      def pvar = it[2] 
      def psam = it[3]
      def traceFile = file("${params.outdir}/channel-trace/unified_plink2_${chr}.txt")
      traceFile.text = "CHR: ${chr}\nPGEN: ${pgen}\nPVAR: ${pvar}\nPSAM: ${psam}\n"
    }
    return it
  }
  .set { unified_plink2_files }

  // Extract pvar files for the existing pipeline
  unified_plink2_files
  .map { chrid, pgen, pvar, psam -> tuple(chrid, pvar) }
  .view { params.dev ? "DEBUG: genotypes_pvar_0: $it" : null }
  .map { 
    if (params.dev) {
      file("${params.outdir}/channel-trace").mkdirs()
      def chr = it[0]
      def pvar = it[1]
      def traceFile = file("${params.outdir}/channel-trace/genotypes_pvar_0_${chr}.txt")
      traceFile.text = "CHR: ${chr}\nPVAR: ${pvar}\n"
    }
    return it
  }
  .set { genotypes_pvar_0 }

  make_geno_pvar_snpid_unique_pvar(genotypes_pvar_0)
  make_geno_pvar_snpid_unique_pvar.out
  .view { params.dev ? "DEBUG: genotypes_pvar_after_unique: $it" : null }
  .map { 
    if (params.dev) {
      file("${params.outdir}/channel-trace").mkdirs()
      def chr = it[0]
      def genFile = it[1]
      def traceFile = file("${params.outdir}/channel-trace/genotypes_pvar_${chr}.txt")
      traceFile.text = "CHR: ${chr}\nFILE: ${genFile}\n"
    }
    return it
  }
  .set { genotypes_pvar }

  // SNPlist
  if (params.snplist) {
    if (params.dev) { println("DEBUG: Using user-provided SNPlist: ${params.snplist}") }
    Channel.fromPath("${params.snplist}", type: 'file')
    .view { params.dev ? "DEBUG: user_snplist_input: $it" : null }
    .set { ch_input_snplist_0 }
    sort_user_snplist(ch_input_snplist_0)
    sort_user_snplist.out
    .view { params.dev ? "DEBUG: sorted_user_snplist: $it" : null }
    .map { 
      if (params.dev) {
        file("${params.outdir}/channel-trace").mkdirs()
        def traceFile = file("${params.outdir}/channel-trace/user_snplist.txt")
        traceFile.text = "SNPLIST_FILE: ${it}\n"
      }
      return it
    }
    .set { ch_input_snplist }
  } else {
    if (params.dev) { println("DEBUG: Creating SNPlist from pvar files") }
    genotypes_pvar.map {x,y -> y}.collect()
    .view { params.dev ? "DEBUG: pvar_files_for_snplist: $it" : null }
    .map { 
      if (params.dev) {
        file("${params.outdir}/channel-trace").mkdirs()
        def traceFile = file("${params.outdir}/channel-trace/pvar_files_for_snplist.txt")
        traceFile.text = "PVAR_FILES: ${it.join('\n')}\n"
      }
      return it
    }
    .set { pvar_files_collected }
    make_snplist_from_pvar(pvar_files_collected)
    make_snplist_from_pvar.out
    .view { params.dev ? "DEBUG: generated_snplist: $it" : null }
    .map { 
      if (params.dev) {
        file("${params.outdir}/channel-trace").mkdirs()
        def traceFile = file("${params.outdir}/channel-trace/generated_snplist.txt")
        traceFile.text = "GENERATED_SNPLIST: ${it}\n"
      }
      return it
    }
    .set { ch_input_snplist }
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
  .view { params.dev ? "DEBUG: filter_bad_values_2_output: $it" : null }
  .map { 
    if (params.dev) {
      file("${params.outdir}/channel-trace").mkdirs()
      def chr = it[0]
      def sumstatFile = it[1]
      def traceFile = file("${params.outdir}/channel-trace/filter_bad_values_2_${chr}.txt")
      traceFile.text = "CHR: ${chr}\nFILE: ${sumstatFile}\n"
    }
    return it
  }
  .join(genotypes_pvar)
  .view { params.dev ? "DEBUG: after_join_genotypes_pvar: $it" : null }
  .combine(ch_input_snplist)
  .view { params.dev ? "DEBUG: after_combine_snplist: $it" : null }
  .join(ch_ldfiles)
  .view { params.dev ? "DEBUG: ch_to_map_final: $it" : null }
  .map { 
    if (params.dev) {
      file("${params.outdir}/channel-trace").mkdirs()
      def chr = it[0]
      def file1 = it[1]
      def file2 = it[2]  
      def snplist = it[3]
      def ldfile1 = it[4]
      def ldfile2 = it[5]
      def traceFile = file("${params.outdir}/channel-trace/ch_to_map_${chr}.txt")
      traceFile.text = "CHR: ${chr}\nSUMSTAT: ${file1}\nPVAR: ${file2}\nSNPLIST: ${snplist}\nLD_BIN: ${ldfile1}\nLD_INFO: ${ldfile2}\n"
    }
    return it
  }
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
  unified_plink2_files
}

workflow wf_sbayesr_calc_score {

  take:
  ch_formatted_posteriors
  variant_maps_for_sbayesr
  sumstat
  unified_plink2_files

  main:

  // Use the unified PLINK2 files for scoring
  unified_plink2_files
  .view { params.dev ? "DEBUG: genotypes_for_scoring: $it" : null }
  .set { genotypes_0 }

  make_geno_pvar_snpid_unique_pvar_psam_pgen(genotypes_0)
  make_geno_pvar_snpid_unique_pvar_psam_pgen.out.set { genotypes }

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
 // if(params.concat_genotypes){

 //   // concat all genotypes
 //   genotypes
 //   .map { chr, bed, bim, fam -> 
 //       def canonicalBed = new File(bed.toString()).canonicalPath
 //       def canonicalBim = new File(bim.toString()).canonicalPath
 //       def canonicalFam = new File(fam.toString()).canonicalPath
 //       return "$canonicalBed $canonicalBim $canonicalFam"
 //   }
 //   .collectFile(){ content -> [ "allgenotypes.txt", content + '\n' ] }
 //   .set { ch_all_genofile }
 //   concat_genotypes(ch_all_genofile)

 //   // join posteriors and genotypes all chromosomes
 //   ch_concatenated_posteriors
 //   .join(concat_genotypes.out)
 //   .set{ ch_calc_score_input_all }
 // }else{
 //   Channel.empty().set { ch_calc_score_input_all }
 // }

  Channel.empty().set { ch_calc_score_input_all }

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

