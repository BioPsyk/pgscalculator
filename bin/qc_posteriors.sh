#!/bin/bash

rscript="${1}"
posteriors="${2}"
inputfile="${3}"
outprefix="${4}"

Rscript ${rscript} ${posteriors} ${inputfile} ${outprefix}


