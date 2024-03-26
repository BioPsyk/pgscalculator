#!/bin/bash

input1="$1"
input2="$2"
outfile="$3"

paste <(zcat ${input2}) <(zcat ${input1}) | cut -f1-2,4- | gzip -c >  ${outfile}

