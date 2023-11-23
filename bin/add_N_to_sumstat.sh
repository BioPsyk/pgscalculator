#!/bin/bash

file="$1"
meta="$2"
priority="$3"

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

function add_totalN_from_Nca_Nco(){
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
  NR > 1 {
    
    # check if chr is NA (which can happen after reverse build liftover)
    if ( $1 == "NA" ){
        print FNR, "chr is NA" > "skipped_rows"
        next
    }

    caseN = $(caseIdx)
    controlN = $(controlIdx)
  
    # Check if caseN and controlN are numbers
    if (caseN !~ /^[0-9]+$/ || controlN !~ /^[0-9]+$/) {
        print FNR, "Non-numeric value in row" > "skipped_rows"
        next
    }
  
    ntot = caseN + controlN
  
    if (existingIdx) {
        $(existingIdx) = int(ntot)
        print $0
    } else {
        print $0, int(ntot)
    }
  }
  ' <(zcat ${file})
elif $tfcol_N ;then
  awk -vOFS="\t" '
  NR==1{
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^N$/) {
	existingIdx=i
      }
    }
    print $0; next
  }
  NR>1{
    $(existingIdx)=int($(existingIdx))
    print $0
  }
  ' <(zcat ${file})

}

function add_effectiveN_from_Nca_Nco(){
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
  NR > 1 {
    
    # check if chr is NA (which can happen after reverse build liftover)
    if ( $1 == "NA" ){
        print FNR, "chr is NA" > "skipped_rows"
        next
    }

    caseN = $(caseIdx)
    controlN = $(controlIdx)
  
    # Check if caseN and controlN are numbers
    if (caseN !~ /^[0-9]+$/ || controlN !~ /^[0-9]+$/) {
        print FNR, "Non-numeric value in row" > "skipped_rows"
        next
    }
  
    neff = 4 / ((1 / caseN) + (1 / controlN))
  
    if (existingIdx) {
        $(existingIdx) = int(neff)
        print $0
    } else {
        print $0, int(neff)
    }
  }
  ' <(zcat ${file})
elif $tfcol_N ;then
  awk -vOFS="\t" '
  NR==1{
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^N$/) {
	existingIdx=i
      }
    }
    print $0; next
  }
  NR>1{
    $(existingIdx)=int($(existingIdx))
    print $0
  }
  ' <(zcat ${file})

}

function use_a_stats_value_and_add_as_N_column(){
  awk -vOFS="\t" -vneff="${1}" 'NR==1{print $0,"N"}; NR>1{print $0, int(neff)}' <(zcat ${file})
}

#priority="effectiveN"
#priority="totalN"
if ${priority}=="totalN"; then
  if ${tfcol_N};then
    #use_variant_specifific_totalN, which is the one already present
    zcat ${file}
  elif ${tfcol_CaseN} && ${tfcol_ControlN} ;then
    add_totalN_from_Nca_Nco
  elif ${tfTotalN};then
    use_a_stats_value_and_add_as_N_column "${TotalN}"
  else
    echo "must provide the col_N, col_CaseN and col_ControlN, or stats_EffectiveN"
    exit 1
  fi
else if ${priority}=="effectiveN"; then
  if ${tfcol_CaseN} && ${tfcol_ControlN} ;then
    add_effectiveN_from_Nca_Nco
  elif [[ "${tfEffectiveN}" != "" ]];then
    use_a_stats_value_and_add_as_N_column ${EffectiveN}
  else
    echo "must provide the col_N, col_CaseN and col_ControlN, or stats_EffectiveN"
    exit 1
  fi
else
  echo "priority option does not exist: ${priority}"
  exit 1
fi

