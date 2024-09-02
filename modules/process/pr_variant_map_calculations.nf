// Nextflow processes to make initial variant map

process sort_user_snplist {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     

    input:
    path("snplist")

    output:
    path("snplist_sorted")

    script:
        """
        LC_ALL=C sort -u -k1,1 "snplist" > "snplist_sorted"
        """
}

process make_snplist_from_bim {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     
    memory '2 GB'


    input:
    path(bims)

    output:
    path("snplist_sorted")

    script:
        """
        for chrfile in ${bims}
        do
          awk '{print \$2}' \$chrfile >> "snplist"
        done
        LC_ALL=C sort -u -k1,1 "snplist" > "snplist_sorted"
        """
}

process variant_map_for_prscs {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     

    input:
    tuple val(chr), path(ss), path(bim), path(ldmeta)

    output:
    tuple val(chr), path("${chr}_variants_mapfile")

    script:
        """

        ### This is outdated, needs to be updated

        # b37 b38 a1 a2 markername
        # chr:pos chr:pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" 'NR>1{print \$1":"\$2, \$3":"\$4, \$7, \$8, \$6 }' ${ss} > ss2
        # chr pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" '{print \$1":"\$4, \$5, \$6, \$2 }' ${bim} > bim2
        # chr pos a1 a2 markername
        awk -vFS=" " -vOFS="\t" 'NR>1{print \$1":"\$3, \$4, \$5, \$2 }' ${ldmeta} > ld2

        # use same variant map as for sbayesr
        variant_map_for_sbayesr.sh ss2 bim2 ld2 "${chr}_variants_mapfile"
        """
}

process variant_map_for_sbayesr {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     

    input:
    tuple val(chr), path(ss), path(bim), path(snplist), path(ldbin), path(ldinfo)

    output:
    tuple val(chr), path("${chr}_variants_mapfile"), path("${chr}_variants_mapfile_no_NA"), emit: map

    script:
        """

        # b37 b38 a1 a2 markername
        # chr:pos chr:pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" 'NR>1{print \$1":"\$2, \$3":"\$4, \$7, \$8, \$6 }' ${ss} > ss2
        # chr pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" '{print \$1":"\$4, \$5, \$6, \$2 }' ${bim} > bim2
        # chr pos a1 a2 markername
        awk -vFS=" " -vOFS="\t" 'NR>1{print \$1":"\$4, \$5, \$6, \$2 }' ${ldinfo} > ld2
        # markername
        cp ${snplist} snp2

        variant_map_for_sbayesr.sh ss2 bim2 snp2 ld2 "${params.gbuild}" "${params.lbuild}" "${chr}_variants_mapfile"
        awk -vFS="\t" -vOFS="\t" '{if (\$9 != "NA") { print \$0 }}' "${chr}_variants_mapfile"  > "${chr}_variants_mapfile_no_NA"
        """
}

// Concatenate variant_map
process concatenate_variant_map {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true


    input:
        path(chrfiles)
    
    output:
        tuple val("all"), path("variant_map.gz")
    script:
        """
        echo "b37     b38     ss_SNP  ss_A1   ss_A2   bim_SNP bim_A1  bim_A2  ld_SNP  ld_A1   ld_A2"  > "variant_map"
        for chrfile in ${chrfiles}
        do
          tail -n+2 \$chrfile >> "variant_map"
        done
        gzip variant_map
        """
}

process filter_sumstat_variants_on_map_file {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     

    input:
    tuple val(chr), path(ss), path("map"), path("map_noNA")

    output:
    tuple val(chr), path("${chr}_subset_sumstats_map"), path("${chr}_subset_sumstats_map_noNA")

    script:
        """
        # subset on chr:pos:a1:a2
        subset_sumstat_on_mapfile.sh ${ss} ${map} "${chr}_subset_sumstats_map"
        subset_sumstat_on_mapfile.sh ${ss} ${map_noNA} "${chr}_subset_sumstats_map_noNA"
        """
}

