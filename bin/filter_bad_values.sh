#!/bin/bash

file=${1}

# Define output files
exclusion_file="excluded"

# Process the file
awk -vFS="\t" -vOFS="\t" '
{
  if (NR == 1) {
      # Process header row
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

  # Extract necessary fields
  beta = $(b_col);
  se = $(se_col);
  eaf = $(eaf_col);

  # Check exclusion criteria
  if (beta == 0 || se == 0 || eaf == 0 || eaf == 1) {
      print $0 > "'$exclusion_file'";
  } else {
      print $0;
  }
}' <(zcat "${file}")

