#!/bin/bash

# Check number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 sumstat mapfile maffile posterior benchmark"
    exit 1
fi

# Define your file paths
sumstat="$1"
mapfile="$2"
maffile="$3"
posterior="$4"
benchmark="$5"

# Target column names in the order you want them
colnames="RSID CHR POS EffectAllele OtherAllele B SE Z P"

# Process the sumstat file to join information from map, posterior, and benchmark
awk -v FS="\t" -v OFS="\t" \
 -v colnames="$colnames" \
 -v mapfile="$mapfile" \
 -v maffile="$maffile" \
 -v posterior="$posterior" \
 -v benchmark="$benchmark" \
'
    BEGIN {

        # Read mapfile into map array (ss_SNP to bim_SNP)
        while ((getline < mapfile) > 0) map[$3] = $6;

        # Read maffile
        FS=" ";
        while ((getline < maffile) > 0) mafHash[$2] = $5;
        FS="\t";

        # Read posterior and benchmark files into their respective arrays
        while ((getline < posterior) > 0) posteriorHash[$4] = $3;
        while ((getline < benchmark) > 0) benchmarkHash[$4] = $3;
        
        # Determine the index of each target column
        n = split(colnames, cols, " ");
        for (i=1; i<=n; i++){
          colIdx[cols[i]] = -1;
          printf "%s%s", cols[i], OFS ;
        }
        print "genoID" OFS "MAF" OFS "postEffect" OFS "benchEffect"
        for (i=1; i<=NF; i++) {
          if ($i in colIdx){
            colIdx[$i] = i;
          }
        }
        getline
        for (i=1; i<=NF; i++) {
          if ($i in colIdx){
            colIdx[$i] = i;
          }
        }
    }
    FNR > 1 { 
        # Assuming RSID will always be in column 2
        rsid = $4; 
        if (rsid in map) {
            genoID = map[rsid];
        } else {
            print "ssID not in map: " rsid > "/dev/stderr";
            exit 1;
        }

        postEffect = (genoID in posteriorHash) ? posteriorHash[genoID] : "NA";
        benchEffect = (genoID in benchmarkHash) ? benchmarkHash[genoID] : "NA";
        mafVal = (genoID in mafHash) ? mafHash[genoID] : "NA";
        
        # Print the selected columns for this row
        for (i=1; i<=n; i++) {
            idx = colIdx[cols[i]];
            printf "%s%s", (idx > 0 ? $idx : "NA"), (i<n ? OFS : "");
        }
        print OFS genoID OFS mafVal OFS postEffect OFS benchEffect; # Append the new fields
    }
' "$sumstat" 

