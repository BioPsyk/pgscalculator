#!/bin/bash

rscript="${1}"
posteriors="${2}"
inputfile="${3}"
inputfile_header="${4}"

Rscript ${rscript} ${posteriors} ${inputfile} ${inputfile_header}

