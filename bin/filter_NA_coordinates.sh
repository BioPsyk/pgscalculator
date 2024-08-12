#!/bin/bash

file="${1}"

# Process the file
awk -vFS="\t" -vOFS="\t" '
$1!="NA" && $2!="NA"{
  print
}' ${file}


