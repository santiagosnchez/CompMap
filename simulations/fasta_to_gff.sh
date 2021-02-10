awk '{ if ($0 ~ /^>/){ col1=substr($0,2) } else { print col1"\t.\ttranscript\t1\t"(length($0))"\t.\t.\t.\tgene_id="(substr(col1,4)) }}' $1
