#!/bin/bash

file="$1"
meta="$2"

#helpers
function selRightHand(){
  echo "${1#*:}"
}
function selColRow(){
  grep ${1} ${2}
}

#recode as true or false
function recode_to_tf(){
  if [ "$1" == "" ]; then
    echo false
  else
    echo true
  fi
}

# What is colname according to meta data file
EAF="$(selRightHand "$(selColRow "^cleansumstats_col_EAF:" $meta)" | sed 's/\s\+//g')"

# Not forwarded right now
#EAF_1KG="$(selRightHand "$(selColRow "^cleansumstats_col_EAF_1KG:" $meta)" | sed 's/\s\+//g')"

# true or false (exists or not)
tfEAF="$(recode_to_tf $EAF)"

# Get header
HEADER=$(zcat "$file" | head -n 1)

if ${tfEAF};then
  # Extra check that EAF exists
  if ! grep -q -w "EAF" <<< "$HEADER"; then
    echo "EAF is not present in header despite being present in the metafile"
    exit 1
  fi
  zcat $file
else 
  # Check for EAF_1KG in header and change it to EAF if exists
  if [[ "$HEADER" == *"EAF_1KG"* ]]; then
    zcat "$file" | sed '1 s/EAF_1KG/EAF/'
  else
    echo "Neither EAF nor EAF_1KG in header of sumstat"
    exit 1
  fi
fi

