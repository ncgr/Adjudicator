#!/bin/bash

if [ $# -lt 4 ];then
    echo "USAGE:./$0 <PREFIX> <GFF3> <GFF3> <gfa.tsv>"
    exit 1
fi

name=$1
input1=$2
input2=$3
gfa=$4

bedtools intersect -wao -a <( awk -F"\t" '{if($3 == "gene"){print}}' $input1 ) -b <( awk -F"\t" '{if($3 == "gene"){print}}' $input2 ) > ${name}.wao.gff3

bedtools intersect -v -a <( awk -F"\t" '{if($3 == "gene"){print}}' $input2 ) -b <( awk -F"\t" '{if($3 == "gene"){print}}' $input1 ) > ${name}.unique.gff3

python3 collapse_overlapping.py ${name}.wao.gff3 ${gfa} ${gfa} 0.00000001 0 ${name} | sort -u | awk -F"\t" '{if(a[$1":"$4":"$5]){}else{a[$1":"$4":"$5]=1;print}}' > ${name}.collapsed.gff3

cat ${name}.collapsed.gff3 ${name}.unique.gff3 | sort -u | awk -F"\t" '{if(a[$1":"$4":"$5]){}else{a[$1":"$4":"$5]=1;print}}' | sort -k1,1 -k4,4n > ${name}.collapsed.final.gff3
