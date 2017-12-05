#!/bin/bash

#  gather_counts_stdin_var.sh
#
#   takes a list of files from stdin
#   returns a table of count values to stdout
#   automatically determines the size of the source files 
#
#   usage: find -maxdepth 4 -name *gencode*body* | 
#            gather_counts_stdin_6col.sh > gencode_genebody_chipseq_counts.txt
#                
#
# James Purzner
# Stanford University, Department of Developmental Biology 
# Fuller Lab 
# Oct 24, 2016


touch all_counts.txt

DATA=$(cat) 




all_cols="$(echo  "${DATA[@]}" | xargs head -1 | grep -P '\t' | awk '{print NF}' | uniq -c | sort -r | head -1 | xargs | cut -d" " -f2-)"

(>&2 echo "autodetected $all_cols columns")

n=1

while read f 
do
    # added to handle empty or incomplete files
    cols=$(head -n 1 $f | awk '{print NF}')
    rows=$(wc -l $f | awk '{print $1}')
    (>&2 echo $f)
    (>&2 echo $cols)
    (>&2 echo $rows)
    if [[ $rows -gt 0 ]] && [[ $cols -eq $all_cols ]]; then
	    if [[ $n -eq 1 ]]; then
		case "$cols" in
		1) printf "gene\t" >> all_counts.txt
		;;
		2) printf "name1\tname2\t" >> all_counts.txt
		;;
		3) printf "chr\tstart\tstop\t" >> all_counts.txt
		;;
		4) printf "chr\tstart\tstop\tname\t" >> all_counts.txt
		;;
		5) printf "chr\tstart\tstop\tname\t" >> all_counts.txt
		;;
		6) printf "chr\tstart\tstop\tname\tempty\t" >> all_counts.txt
		;;
		7) printf "chr\tstart\tstop\tname\tempty\tstrand\t" >> all_counts.txt
		;;
		esac
		colname=${f}
		colname="${colname:2:${#colname}-1}"
		printf "%s\n" "${colname}" >> all_counts.txt
		(>&2 echo "$colname")
		cat $f >> all_counts.txt
	    else
		touch tmp.txt
		touch tmp2.txt
		colname=${f}
		colname="${colname:2:${#colname}-1}"
		printf "%s\n" "${colname}" >> tmp.txt
		(>&2 echo "$colname")
		cat $f | awk -v var=$all_cols 'BEGIN {OFS="\t"; print $var}' >> tmp.txt
		paste all_counts.txt tmp.txt >> tmp2.txt
		mv -f tmp2.txt all_counts.txt
		rm tmp.txt
		(>&2 echo $n)
	    fi
	    n=$(( $n + 1 ))
    fi 

done <<< "$DATA"

cat all_counts.txt
rm all_counts.txt

# misc code 
#colname="${colname:2:${#colname}-1}"
