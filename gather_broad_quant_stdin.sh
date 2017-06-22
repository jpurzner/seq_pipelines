#!/bin/bash

#  gather_counts_find.sh
#
#   takes a list of files from stdin
#   returns a table of count values to stdout
#
#   usage: find . -maxdepth 2 -name RNAseq_counts_mm9_ucsc_canonical_12.bed |
#                bash gather_counts_find.sh > out.txt
#
#


touch all_counts.txt

n=1
while read f 
do
    
    cols=$(head -n 1 $f | awk '{print NF}')
    rows=$(wc -l $f | awk '{print $1}')
    (>&2 echo $f)
    (>&2 echo $cols)
    (>&2 echo $rows)
    if [[ $rows -gt 0 ]] && [[ $cols -eq 8 ]]; then
	    if [[ $n -eq 1 ]]; then
		printf "chr\tstart\tstop\tname\tempty\tstrand\t" >> all_counts.txt
		colname=${f}
		colname=$(basename ${colname})
		colname=${colname%_broad_quant.bed}
		(>&2 echo "$colname")
		printf "%s\t%s\n" "${colname}" "${colname}.input"  >> all_counts.txt
		cat $f | awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8}' >> all_counts.txt
	    else
		touch tmp.txt
		touch tmp2.txt
		colname=${f}
		colname=$(basename ${colname})
		colname=${colname%_broad_quant.bed}
		(>&2 echo "$colname")
		printf "%s\t%s\n" "${colname}" "${colname}.input" >> tmp.txt
		cat $f | awk '{OFS="\t"; print $7, $8}' >> tmp.txt
		paste all_counts.txt tmp.txt >> tmp2.txt
		mv -f tmp2.txt all_counts.txt
		rm tmp.txt
		(>&2 echo $n)
	    fi
	    n=$(( $n + 1 ))
    fi 

done

cat all_counts.txt
rm all_counts.txt

# misc code 
#colname="${colname:2:${#colname}-1}"
