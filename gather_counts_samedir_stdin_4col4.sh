#!/bin/bash

#  gather_counts_find.sh
#
#   takes a list of files from stdin
#   returns a table of count values to stdout
#
#   usage: find *_mm9_ucsc_canonical_12.bed |
#                bash gather_counts_find.sh > out.txt
#
#




touch all_counts.txt

n=1
while read f 
do
    if [ $n -eq 1 ]; then
        printf "name\tchr\tstart\tstop\t" >> all_counts.txt
        colname=${f}
	printf "%s\n" "${colname}" >> all_counts.txt
        cat $f | awk '{OFS="\t"; print $1, $2, $3,  $4, $5}' >> all_counts.txt
    else
        touch tmp.txt
        touch tmp2.txt
	colname=${f}
        printf "%s\n" "${colname}" >> tmp.txt
        cat $f | awk '{OFS="\t"; print $5}' >> tmp.txt
        paste all_counts.txt tmp.txt >> tmp2.txt
        mv -f tmp2.txt all_counts.txt
        rm tmp.txt
    fi
    n=$(( $n + 1 ))
done

cat all_counts.txt
rm all_counts.txt