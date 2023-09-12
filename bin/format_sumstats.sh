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

awk -v format="$format" -vFS="\t" -vOFS="\t"'
# Process the map file
NR==FNR {
    if ($1 == format) {
        for (i=2; i<=NF; i++) {
            if ($i != "NA") {
                header[i] = $i
            }
        }
    } else if ($1 == format"_cix") {
        for (i=2; i<=NF; i++) {
            if ($i != "NA") {
                col_idx[header[i]] = $i
            }
        }
    }
    next
}

# Process the input file
FNR==1 {
    for (i=1; i<=NF; i++) {
        input_col[$i] = i
    }
    for (h in header) {
        if (header[h] == "CHR") {
            chr_col = input_col[col_idx[header[h]]]
        }
    }
    next
}

{
    if (chr_col) {
        outfile = "'"$output_prefix"'" "_" $chr_col ".tsv"
    } else {
        outfile = "'"$output_prefix"'.tsv"
    }
    for (h in header) {
        printf $input_col[col_idx[header[h]]] (h < NF ? OFS : "\n") > outfile
    }
}
' "$map_file" "$input_file"
