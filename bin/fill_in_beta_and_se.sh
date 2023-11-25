#!/bin/bash

inputfile=$1

# Process the file
awk -vFS='\t' -vOFS='\t' '
{
    if (NR == 1) {
        # Process header row
        for (i = 1; i <= NF; i++) {
            if ($i == "Z") z_col = i;
            if ($i == "N") n_col = i;
            if ($i == "EAF") eaf_col = i;
            if ($i == "B") b_col = i;
            if ($i == "SE") se_col = i;
        }

        # Print header and add B and SE columns if not present
        print $0, (b_col ? "" : "B"), (se_col ? "" : "SE");
        next;
    }

    # Extract necessary fields using the identified column indices
    z = $(z_col);
    N = $(n_col);
    freq_l = $(eaf_col);

    # Check if B and SE are missing and calculate them
    if (!b_col || !se_col) {
        beta = z / (sqrt(2 * freq_l * (1 - freq_l) * (N + z^2)));
        se = 1 / (sqrt(2 * freq_l * (1 - freq_l) * (N + z^2)));
    }

    # Print the original line with B and SE appended if they were missing
    print $0, (!b_col ? beta : ""), (!se_col ? se : "");
}' ${inputfile}

