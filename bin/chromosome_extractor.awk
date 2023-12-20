# A script to pck out the 1000 first rows for each chromsomes for a sumstat
# This is to produce representative example data

BEGIN { header = "" }

# Capture and print the header
NR == 1 {
    header = $0
    print $0
    next
}

# For each line, count the occurrences of the first field (chromosome)
{
    if (count[$1] < 1000) {
        print $0
        count[$1]++
    }
}

