inbim=$1
outbim=$2

# Create a temporary file for duplicates
awk '{print $3}' $inbim | grep -v '^#' | sort | uniq -d > duplicates.txt

if [ -s duplicates.txt ]  # Check if duplicates file has content
then
  awk -v OFS="\t" '
    FNR==NR {  # First pass: read duplicates file
      dup[$1]=1;
      next;
    }
    {
      if ($1 ~ /^#/) {
        print;
        next;
      }
      if ($3 in dup) {
        count[$3]++;
        $3 = $3 "_" count[$3];
      }
      print;
    }' duplicates.txt $inbim > $outbim
else
  cp $inbim $outbim
fi
