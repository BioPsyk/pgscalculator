#!/bin/bash

input_file="${1}"
output_prefix="${2}"

awk -voutput_prefix="${output_prefix}" -vFS="\t" -vOFS="\t" \
'
BEGIN {
  # Read the first line of the main input file
  getline

  # Store infile header name positions
  for (i=1; i<=NF; i++) {
      input_col[$i] = i

      # Check if "CHR" is available in the header
      if ($i == "CHR") {
        chr_col=i
      }
  }
  header=$0

  # init header in all output files, Loop over chromosome numbers 1-22
  for (chr = 1; chr <= 22; chr++) {
    outfile = sprintf("%s_%s.tsv", output_prefix, chr)
    print header > outfile
  }
}
{
    print $0 > output_prefix "_" $(chr_col) ".tsv"
}
' <(zcat ${input_file})

