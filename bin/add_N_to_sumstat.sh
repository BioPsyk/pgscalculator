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

#what is colname according to meta data file
TraitType="$(selRightHand "$(selColRow "^stats_TraitType:" $meta)" | sed 's/\s\+//g')"
TotalN="$(selRightHand "$(selColRow "^stats_TotalN:" $meta)" | sed 's/\s\+//g')"
CaseN="$(selRightHand "$(selColRow "^stats_CaseN:" $meta)" | sed 's/\s\+//g')"
ControlN="$(selRightHand "$(selColRow "^stats_ControlN:" $meta)" | sed 's/\s\+//g')"
EffectiveN="$(selRightHand "$(selColRow "^stats_EffectiveN:" $meta)" | sed 's/\s\+//g')"
col_N="$(selRightHand "$(selColRow "^cleansumstats_col_N:" $meta)" | sed 's/\s\+//g')"
col_CaseN="$(selRightHand "$(selColRow "^cleansumstats_col_CaseN:" $meta)" | sed 's/\s\+//g')"
col_ControlN="$(selRightHand "$(selColRow "^cleansumstats_col_ControlN:" $meta)" | sed 's/\s\+//g')"

#true or false (exists or not)
tfTraitType="$(recode_to_tf $TraitType)"
tfTotalN="$(recode_to_tf $TotalN)"
tfCaseN="$(recode_to_tf $CaseN)"
tfControlN="$(recode_to_tf $ControlN)"
tfEffectiveN="$(recode_to_tf $EffectiveN)"
tfcol_N="$(recode_to_tf $col_N)"
tfcol_CaseN="$(recode_to_tf $col_CaseN)"
tfcol_ControlN="$(recode_to_tf $col_ControlN)"


if ${tfcol_CaseN} && ${tfcol_ControlN} ;then
  
  awk -vOFS="\t" '
  NR==1 {
    
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^CaseN$/) {  
        caseIdx = i
      } else if ($i ~ /^ControlN$/) {  
        controlIdx = i
      } else if ($i ~ /^N$/) {
        existingIdx = i
      }
    
    }
    if (!caseIdx || !controlIdx) {
      print "Error: Required columns caseN or controlN not found in the header." > "/dev/stderr"
      exit 1
    }
    if (existingIdx) {
      print $0 
    }else {
      print $0, "N"
    }
  }
  NR>1 {
    caseN = $(caseIdx)
    controlN = $(controlIdx)
    neff = 4 / ((1 / caseN) + (1 / controlN))
    if (existingIdx) {
      $(existingIdx)=neff
    }
    print $0
  }
  ' <(zcat ${file})
elif $tfcol_N ;then
  zcat ${file}
elif [[ "${tfEffectiveN}" != "" ]];then
  awk -vOFS="\t" -vneff="${EffectiveN}" 'NR==1{print $0,"N"}; NR>1{print $0, neff}' <(zcat ${file})
else
  echo "must provide the col_N, col_CaseN and col_ControlN, or stats_EffectiveN"
  exit 1
fi

