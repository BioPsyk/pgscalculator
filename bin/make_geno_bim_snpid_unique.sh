inbim=$1
outbim=$2

awk '{print $2}' $inbim | sort | uniq -d > duplicates.txt

if [ "$(head -n1 duplicates.txt | wc -l)" == "1" ]
then
  awk 'BEGIN {
      FS=OFS="\t";
      # Load potential duplicates to identify them quickly
      while (getline < "duplicates.txt" > 0)
          dups[$1] = 1;
  }
  {
      if ($2 in dups) {
          id = $1 ":" $4 ":" $5 ":" $6;  # Construct the ID based on chr:pos:a1:a2
          count[id]++;  # Increment the counter for this ID
          if (count[id] > 1) {  # If this ID has appeared before, append a counter
              $2 = id ":" count[id];
          } else {
              $2 = id;  # First occurrence of this ID, use without appending a number
          }
      }
      print;  # Print the current line with possibly modified ID
  }' $inbim > $outbim

  awk '{print $2}' $outbim | sort | uniq -d > duplicates2.txt
  if [ "$(head -n1 duplicates2.txt | wc -l)" == "1" ]
  then
    echo "there are still duplicated genotypes in the bim file"
    exit 1
  else
    :
  fi
else
  ln -s  $inbim $outbim
fi
