#!/bin/bash

input="$1"
torm="$2"
outfile="$3"

if [ "${torm}" == "2" ]; then
  cut -f1-2,5- ${input} > ${outfile}
elif [ "${torm}" == "1" ]; then
  cut -f3- ${input} > ${outfile}
else
  echo "This torm option does not exist."
fi

