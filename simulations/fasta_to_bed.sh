awk '{ if ($0 ~ /^>/){ col1=substr($0,2) } else { print col1"\t"0"\t"(length($0)) }}' $1
