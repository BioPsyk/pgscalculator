#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file> <map_file> <format>"
    exit 1
fi

input_file="$1"
map_file="$2"
format="$3"

# Methods specific delim
MOFS="\t"
if [ ${format} == "sbayesr" ]; then
  MOFS=" "
fi

awk -v format="$format" -v mapfile="$map_file" -vFS="\t" -vOFS="${MOFS}" \
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

  # Print the header based on ordered_cols and out_header
  for (i in ordered_cols) {
      printf out_header[ordered_cols[i]] (i < length(ordered_cols) ? OFS : "\n")
  }

}

{
  for (i = 1; i < length(ordered_cols); i++) {
    printf $(input_col[in_header[ordered_cols[i]]]) OFS 
  }
  printf $(input_col[in_header[ordered_cols[i]]]) "\n"
}
' ${input_file}

