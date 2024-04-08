#!/bin/bash

# Usage: ./find_col_indices.sh "file" "column1,column2,column3"

file="$1"
cols="$2"

# Check if cols are integers (indices) or strings (names)
if [[ $cols =~ ^[0-9,]+$ ]]; then
    # If cols are integers, just print them
    echo $cols
else
    # If cols are names, find and print their indices
    awk -v cols="$cols" -v FS="\t" '
    BEGIN {
        n = split(cols, names, ",");
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            for (j = 1; j <= n; j++) {
                if ($i == names[j]) {
                    idx[j] = i;  # Store index by order
                    found[j] = 1;  # Mark as found
                }
            }
        }
        for (j = 1; j <= n; j++) {
            if (!found[j]) {
                print "Error: Column " names[j] " not found" > "/dev/stderr";
                exit 1;
            }
            printf "%d", idx[j];
            if (j < n) printf ",";
        }
        printf "%s", "\n";
        exit;
    }' "$file"
fi

