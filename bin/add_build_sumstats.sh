#!/bin/bash

input1="$1"
input2="$2"
outfile="$3"

#paste <(zcat ${input2}) <(zcat ${input1}) | cut -f1-2,4- | gzip -c >  ${outfile}

# temporary fix until cleansumstats fixed the empty values for NAs in the GRCh37 map
paste <(zcat ${input2}) <(zcat ${input1}) | awk -vFS="\t" -vOFS="\t" '{for(i=1; i<=NF; i++) if($i=="") $i="NA"; print}' | cut -f1-2,4- | gzip -c >  ${outfile}


