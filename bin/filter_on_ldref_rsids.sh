#!/bin/bash
sfile=$1
rsids=$2

zcat ${sfile} | awk -vrsids_file="${rsids}" -vFS="\t" -vOFS="\t" ' 
BEGIN{
  while (getline < rsids_file) { 
    rsids[$1];
  }
}
FNR == 1 {
  # Process header row to find the RSID column
  for (i = 1; i <= NF; i++) {
    if ($i == "RSID") {
      rsid_col = i;
      print $0;
      break;
    }
  }
}
FNR > 1 {
  if ($(rsid_col) in rsids) {
    print $0;
  }
}
'

