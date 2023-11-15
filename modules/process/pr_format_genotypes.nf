#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process add_rsid_to_genotypes {
    //publishDir "${params.outdir}/intermediates/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'mod_mem'
    
    input:
        tuple val(chr), path("geno.bed"), path("genoX.bim"), path("geno.fam"), path(rsid_ref)

    output:
        tuple val(chr), path("geno2.pgen"), path("geno2.pvar"), path("geno2.psam")
    
    script:
        """
        
        # add rsid and make list of tokeep
        add_rsid_to_genotypes.sh "genoX.bim" ${rsid_ref} > "geno.bim"

        # only keep variants with rsid, remove dups and indels
        plink2 --bfile geno --make-pgen --extract tokeep --out geno2
        """
}

process concat_genotypes {
    publishDir "${params.outdir}/intermediates/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'mod_mem'
    
    input:
       path("allfiles.txt")

    output:
        tuple val("all"), path("allgeno.pgen"), path("allgeno.pvar"), path("allgeno.psam")
    
    script:
        """
        plink2 --pmerge-list allfiles.txt --make-pgen --multiallelics-already-joined --out allgeno
        """

}
//        oprefix=1
//        while read -r c1 c2 c3; do
//          echo "\$(realpath \$c1)"
//          echo "\$(realpath \$c2)"
//          echo "\$(realpath \$c3)"
//          awk '{seen[\$2]++}; END{for (k in seen){if(seen[k]==1){print k}}}' "\$(realpath \$c2)" > include_file 
//          plink --bed "\$(realpath \$c1)" --bim "\$(realpath \$c2)" --fam "\$(realpath \$c3)" --extract include_file --make-bed --out "prefix_\${oprefix}"
//          echo "prefix_\${oprefix}" >> allfiles_converted.txt
//          ((oprefix++))
//        done < allfiles.txt
//        plink --merge-list allfiles_converted.txt --make-bed --out allgeno

        //oprefix=1
        //while read -r c1 c2 c3; do
        //  echo "\$(realpath \$c1)"
        //  echo "\$(realpath \$c2)"
        //  echo "\$(realpath \$c3)"
        //  prefix=\${c1%.bed}
        //  plink --bfile \${prefix} --recode vcf --out "prefix_\${oprefix}"
        //  echo "prefix_\${oprefix}" >> allfiles_converted.txt
        //  ((oprefix++))
        //done < allfiles.txt

        //exit 1
        //plink2 --pmerge-list allfiles.txt bfile --make-bed --out allgeno
          //plink2 --bfile \${prefix} --make-pgen --out "prefix_\${oprefix}"
        //plink --bfile --merge-list allfiles.txt --make-bed --out allgeno

