#!/bin/bash
sfile=$1
rsids=$2

# Read RSIDs into an array in awk
awk -vFS="\t" -vOFS="\t"' 
NR == FNR {
    rsids[$1];
    next;
}
NR != FNR && FNR == 1 {
  # Process header row to find the RSID column
  for (i = 1; i <= NF; i++) {
    if ($i == "RSID") {
      rsid_col = i;
      print $0;
      break;
    }
  }
}
NR != FNR && FNR > 1 {
  if ($(rsid_col) in rsids) {
    print $0;
  }
}' "$rsids" <(zcat ${sfile})

