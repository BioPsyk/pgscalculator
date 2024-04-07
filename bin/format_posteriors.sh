#!/bin/bash

# Assign command line arguments to variables for clarity
post="$1"       # Path to the posteriors file
postCols="$2"   # Columns to select from the posteriors file (e.g., "2 5 8")
map="$3"        # Path to the map file
mapCols="$4"    # Columns to map from and to (e.g., "3 6")

# Use awk to process the files
awk -vFS="\t" -vOFS="\t" -v postCols="$postCols" -v mapCols="$mapCols" '
BEGIN {
    # Split column specifications into arrays
    n = split(postCols, postArr, ",");
    m = split(mapCols, mapArr, ",");

    # Assign split values to meaningful variable names
    idCol = postArr[1];
    a1Col = postArr[2];
    effectCol = postArr[3];
    mapFromCol = mapArr[1];
    mapToCol = mapArr[2];

    # Print header
    print "posteriorID", "EA", "Effect", "genoID";
}

# Process the map file first to populate an associative array
FNR == NR {
    map[$mapFromCol] = $mapToCol;
    next;
}

# Process the posteriors file
FNR > 1 {  # Skip the header line of the posteriors file
    if ($idCol in map) {
        # Print the selected columns and append the mapped genotype from the map array
        print $idCol, $a1Col, $effectCol, map[$idCol];
    } else {
        # If a mapping is not found, print an error message to stderr and exit
        print "Error: No mapping found for SNP ID", $idCol > "/dev/stderr";
        exit 1;
    }
}' "$map" <( awk -vOFS="\t" '{$1=$1; print}' "$post") 


