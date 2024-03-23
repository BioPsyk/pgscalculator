sumstat=$1
mapfile=$2
outfile=$3

head -n1 ${sumstat} > ${outfile}

awk -vFS="\t" -vOFS="\t" '
    NR==FNR{a[$1":"$4":"$5]; next}
    NR!=FNR{if($1":"$2":"$7":"$8 in a){print $0}}
    ' ${mapfile} ${sumstat} >> ${outfile}

