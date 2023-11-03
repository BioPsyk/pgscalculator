#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_file> <map_file> <format> <output_file_prefix>"
    exit 1
fi

input_file="$1"
map_file="$2"
format="$3"
output_prefix="$4"

awk -v format="$format" -v mapfile="$map_file" -v output_prefix="$output_prefix" -vFS="\t" -vOFS="\t" \
'
BEGIN {
  # Process the map file
  getline < mapfile
  
  for (i=2; i<=NF; i++) {
    in_header[i] = $i
  }

  while (getline < mapfile) {
      if ($1 == format) {
          for (i=2; i<=NF; i++) {
              if ($i != "NA") {
                  out_header[i] = $i
              }else{
                  # Set in_header to NA
                  in_header[i] = "NA"
              }
          }
      } else if ($1 == format"_cix") {
          for (i=2; i<=NF; i++) {
              if ($i != "NA") {
                  #col_idx[out_header[i]] = $i
                  #output_name[$i] = out_header[i]
                  ordered_cols[$i] = i
              }
          }
      }
  }

  # Check if format and format_cix are found
  if (!(length(out_header) > 0 && length(ordered_cols) > 0)) {
      print "could not find the right rows in mapfile"
      exit 1
  }

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

  # check that everything in in_header also exists as input_col
  for (k in in_header) {
    if(in_header[k] != "NA"){
      if (!(in_header[k] in input_col)) {
        print "Header not found in input sumstat columns: " in_header[k]
        exit 1
      }
    }
  }

  if (chr_col){
    # init header in all output files, Loop over chromosome numbers 1-22
    for (chr = 1; chr <= 22; chr++) {
      outfile = sprintf("%s_%s.tsv", output_prefix, chr)
      for (i = 1; i < length(ordered_cols); i++) {
        printf out_header[ordered_cols[i]] OFS > outfile
      }
      print out_header[ordered_cols[length(ordered_cols)]] > outfile
    }
  }else{
    ## init header in the only output file, using chr=0
    #chr=0
    #for (h in header) {
    #  outfile = sprintf("%s_%s.tsv", output_prefix, chr)
    #  printf $(input_col[col_idx[header[h]]]) (h < NF ? OFS : "\n") > outfile
    #}
  }
}

{
  if(chr==0){ 
    for (i = 1; i < length(ordered_cols); i++) {
      printf $(input_col[in_header[ordered_cols[i]]]) OFS > output_prefix "_0" ".tsv"
    }
    printf $(input_col[in_header[ordered_cols[i]]]) "\n" > output_prefix "_0" ".tsv"
    
  }else{
    for (i = 1; i < length(ordered_cols); i++) {
      printf $(input_col[in_header[ordered_cols[i]]]) OFS > output_prefix "_" $(chr_col) ".tsv"
    }
    printf $(input_col[in_header[ordered_cols[i]]]) "\n" >  output_prefix "_" $(chr_col) ".tsv"
  }
}
' <(zcat "$input_file")

