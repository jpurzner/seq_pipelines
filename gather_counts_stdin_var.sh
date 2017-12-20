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

#oIFS="$IFS"
#IFS=" "
touch all_counts.txt

DATA=$(cat) 
# need to fix this 
#common=($(common_string.py ${DATA[@]}))



#IFS=$oIFS

all_cols="$(echo  "${DATA[@]}" | xargs head -1 | grep -P '\t' | awk '{print NF}' | uniq -c | sort -r | head -1 | xargs | cut -d" " -f2-)"
(>&2 echo "autodetected $all_cols columns")


# this still needs work use loop that stops when len 1 
common_txt="$(echo  "${DATA[@]}" |  sed -e 'N;s/^\(.*\).*\n\1.*$/\1/')"
common_txt="$(echo  "${common_txt[@]}" |  sed -e 'N;s/^\(.*\).*\n\1.*$/\1/')"
common_txt="$(echo  "${common_txt[@]}" |  sed -e 'N;s/^\(.*\).*\n\1.*$/\1/')"
(>&2 echo  ${#common_txt[@]})


#common_txt="$(echo  "${DATA[@]}" | awk '{print $0}'  | grep -zoP '^(.*)(?=.*?\n\1)')"
(>&2 echo "common text: ${common_txt}")
#(>&2 echo "${common_txt}")




n=1

while read f 
do
    # added to handle empty or incomplete files
    cols=$(head -n 1 $f | awk '{print NF}')
    rows=$(wc -l $f | awk '{print $1}')
    #(>&2 echo $f)
    #(>&2 echo $cols)
    #(>&2 echo $rows)
    if [[ $rows -gt 0 ]] && [[ $cols -eq $all_cols ]]; then
	    if [[ $n -eq 1 ]]; then
		case "$all_cols" in
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
		#colname="${colname:2:${#colname}-1}"
		colname="${colname//$common_txt/}"
		colname="${colname%.*}"
		printf "%s\n" "${colname}" >> all_counts.txt
		(>&2 echo "$colname")
		(>&2 echo "$common_txt")
		cat $f >> all_counts.txt
	    else
		touch tmp.txt
		touch tmp2.txt
		colname=${f}
		#colname="${colname:2:${#colname}-1}"
		colname="${colname//$common_txt/}"
		colname="${colname%.*}"
		printf "%s\n" "${colname}" >> tmp.txt
		(>&2 echo "$colname")
		cat $f | awk -v counts=$all_cols '{OFS="\t"; print $counts}' >> tmp.txt
		paste all_counts.txt tmp.txt >> tmp2.txt
		mv -f tmp2.txt all_counts.txt
		rm tmp.txt
		(>&2 echo $n)
	    fi
	    n=$(( $n + 1 ))
    fi 

done <<< "$DATA"
#cat all_counts.txt | awk 'NR == 1' | sed  's|'$common_txt'||' 
cat all_counts.txt

rm all_counts.txt

