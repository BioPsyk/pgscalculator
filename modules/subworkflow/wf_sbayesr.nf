#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 add_N_effective
 format_sumstats
 force_EAF_to_sumstat
 add_B_and_SE
} from '../process/pr_format_sumstats.nf'
include { 
 calc_posteriors_sbayesr 
 concatenate_sbayes_posteriors
} from '../process/pr_calc_posteriors.nf'
include { calc_score } from '../process/pr_calc_score.nf'
include { 
add_rsid_to_genotypes 
concat_genotypes
} from '../process/pr_format_genotypes.nf'

workflow wf_sbayesr {

    take:
    input

    main:

    if (params.calc_posterior) {

      // metafile for sumstat
      Channel.fromPath("${params.input}/cleaned_metadata.yaml", type: 'file').set { ch_input_metafile }

      // support files from assets
      if (params.mapfile) { mapfile = file(params.mapfile, checkIfExists: true) }


      // format sumstat
      add_N_effective(input, ch_input_metafile, "${params.whichN}")
      force_EAF_to_sumstat(add_N_effective.out, ch_input_metafile)
      add_B_and_SE(force_EAF_to_sumstat.out)
      format_sumstats(add_B_and_SE.out, mapfile, "sbayesr")
      .flatMap { it }
      .map { file ->
        def parts = file.name.split("_")
        [parts[1].replace(".tsv", ""), file]
      }
      .set { sumstats }

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

      sumstats
      .join(ch_ldfiles)
      .set { ch_calc_posteriors }  

      //ch_calc_posteriors
      calc_posteriors_sbayesr(ch_calc_posteriors).set { ch_calculated_posteriors }

    }else{
      input
      .map { file ->
        def chrWithPrefix = file.getBaseName().split("_")[0]
        def chr = chrWithPrefix.replaceAll("chr", "")
        return tuple(chr, file)
      }.set { ch_calculated_posteriors }
    }

    if(params.calc_score){

      // channel of genotypes
      Channel.fromPath("${params.genofile}")
      .splitCsv(sep: '\t', header: false)
      .map { row -> row.collect { it.trim() } } // Trim whitespace from each field
      .map { row -> tuple(row[0], row[1], file("${params.genodir}/${row[2]}")) }
      .filter { type -> type[1] in ['bed', 'bim', 'fam'] }
      .groupTuple()
      .map { chrid, _, files -> [chrid, *files] }
      .set { genotypes }

      // always add rsids based on our dbsnp reference
      if ("${params.gbuild}" == "37") {
           Channel.fromPath("${params.rsid_ref_37}/*")
           .set { ch_rsid_ref }
      } else if ("${params.gbuild}" == "38") {
           Channel.fromPath("${params.rsid_ref_38}/*")
           .set { ch_rsid_ref }
      } else {
          error "Genome build has to be 37 or 38, now it is ${params.gbuild}"
      }

      ch_rsid_ref.map { file ->
        def chrWithPrefix = file.getBaseName().split("_")[0]
        def chr = chrWithPrefix.replaceAll("chr", "")
        return tuple(chr, file)
      }.set {ch_rsid_ref2}

      genotypes
      .join(ch_rsid_ref2)
      .set { ch_add_rsid_to_genotypes }
      add_rsid_to_genotypes(ch_add_rsid_to_genotypes)


      // concat all posteriors
      calc_posteriors_sbayesr.out.map {x,y -> y}.collect().set {ch_collected_posteriors}
      concatenate_sbayes_posteriors(ch_collected_posteriors)
      concatenate_sbayes_posteriors.out.set { ch_concatenated_posteriors }

      // concat all genotypes
      add_rsid_to_genotypes.out
      .map { chr, bed, bim, fam -> 
          def canonicalBed = new File(bed.toString()).canonicalPath
          def canonicalBim = new File(bim.toString()).canonicalPath
          def canonicalFam = new File(fam.toString()).canonicalPath
          return "$canonicalBed $canonicalBim $canonicalFam"
      }
      .collectFile(){ content -> [ "allgenotypes.txt", content + '\n' ] }
      .set { ch_all_genofile }
      concat_genotypes(ch_all_genofile)

      // join posteriors and genotypes (1/2)
      ch_concatenated_posteriors
      .join(concat_genotypes.out)
      .set{ ch_calc_score_input_all }

      // join posteriors and genotypes (2/2)
      ch_calculated_posteriors
      .join(add_rsid_to_genotypes.out)
      .set{ ch_calc_score_input_per_chr }

      // Calc score
      ch_calc_score_input_per_chr
      .mix(ch_calc_score_input_all)
      .set{ ch_calc_score }
      calc_score(ch_calc_score, "${params.calc_posteriors_sbayesr.score_columns}")
    }
}

