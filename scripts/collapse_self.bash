#!/bin/bash

if [ $# -lt 3 ];then
    echo "USAGE:./$0 <PREFIX> <GFF3> <gfa.tsv>"
    exit 1
fi

name=$1
input=$2
gfa=$3

/falafel/ctc/collapse_alfalfa/sw/bedtools2-2.31.1/bin/bedtools intersect -wao -a <( awk -F"\t" '{if($3 == "gene"){print}}' $input ) -b <( awk -F"\t" '{if($3 == "gene"){print}}' $input ) | awk 'BEGIN {FS="\t"} $9 > $18 {print}' > ${name}.self.wao.gff3

/falafel/ctc/collapse_alfalfa/sw/bedtools2-2.31.1/bin/bedtools intersect -wao -a <( awk -F"\t" '{if($3 == "gene"){print}}' $input ) -b <( awk -F"\t" '{if($3 == "gene"){print}}' $input ) | awk -F"\t" '{a[$9]+=1;b[$9]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}END{for(f in a){if(a[f] == 1){print b[f]}}}' > ${name}.self.unique.gff3

python3 collapse_overlapping.py ${name}.self.wao.gff3 ${gfa} ${gfa} 0.00000001 0 $name | sort -u | awk -F"\t" '{if(a[$1":"$4":"$5]){}else{a[$1":"$4":"$5]=1;print}}' > ${name}.self.collapsed.gff3

cat ${name}.self.collapsed.gff3 ${name}.self.unique.gff3 | sort -u | awk -F"\t" '{if(a[$1":"$4":"$5]){}else{a[$1":"$4":"$5]=1;print}}' | sort -k1,1 -k4,4n > ${name}.self.final.gff3
