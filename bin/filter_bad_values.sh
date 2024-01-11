#!/bin/bash

file="${1}"
which="${2}"

# Define output files
exclusion_file="excluded"

# Process the file
awk -vFS="\t" -vOFS="\t" '
{
  if (NR == 1) {
      # Process header row
      b_col = 0
      se_col = 0
      eaf_col = 0
      for (i = 1; i <= NF; i++) {
          if ($i == "B") b_col = i;
          if ($i == "SE") se_col = i;
          if ($i == "EAF") eaf_col = i;
      }

      # Print headers to both output and exclusion files
      print $0;
      print $0 > "'$exclusion_file'";
      next;
  }

  # Check if exclude row
  exclude=0
  if(b_col){ 
    if($(b_col)==0) exclude=1
    if($(b_col)=="NA") exclude=1
  }
  if(se_col){ 
    if($(se_col)==0) exclude=1
    if($(se_col)=="NA") exclude=1
  }
  if(eaf_col){
    if($(eaf_col)==0) exclude=1
    if($(eaf_col)==1) exclude=1
    if($(eaf_col)=="NA") exclude=1
  }

  # Check exclusion criteria
  if (exclude) {
      print $0 > "'$exclusion_file'";
  } else {
      print $0;
  }
}' ${file}


