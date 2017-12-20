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
#IFS=$'\n'
touch all_counts.txt

DATA=$(cat) 
# need to fix this 
#common=($(common_string.py ${DATA[@]}))



(>&2 echo "--------------------------------------------------------------")
(>&2 echo "--------------------------------------------------------------")

all_cols="$(echo  "${DATA[@]}" | xargs head -1 | grep -P '\t' | awk '{print NF}' | uniq -c | sort -r | head -1 | xargs | cut -d" " -f2-)"
(>&2 echo "autodetected $all_cols columns")
file_num=$(echo "$DATA" | wc -l)

(>&2 echo "processing  $file_num bed files")


# identify common string
common_txt=$( printf '%s' "$DATA" |  sed -e 'N;s/^\(.*\).*\n\1.*$/\1/') 
common_txt=$( printf '%s' "$common_txt" |  sed -e 'N;s/^\(.*\).*\n\1.*$/\1/') 
#(>&2 echo  "$common_txt")

com_count=$(echo "$common_txt" | wc -l)
#(>&2 echo "current common length $com_count")

while [[ $com_count -gt 1 ]];
do
    common_txt=$( printf '%s' "$common_txt" |  sed -e 'N;s/^\(.*\).*\n\1.*$/\1/') 
    com_count=$(echo "$common_txt" | wc -l)
    #(>&2 echo "current common length $com_count")
done


#(>&2 echo  "done")
(>&2 echo  "common text: $common_txt")
(>&2 echo "--------------------------------------------------------------")


# create 2 temporary files 
tfile1=$(mktemp tmp1.XXXXXXXXX)
tfile2=$(mktemp tmp2.XXXXXXXXX)
tcount=$(mktemp tcount.XXXXXXXXX)


n=1

while read f 
do
    # added to handle empty or incomplete files
    cols=$(head -n 1 $f | awk '{print NF}')
    rows=$(wc -l $f | awk '{print $1}')
    (>&2 echo "file: $f")
    (>&2 echo "columns: $cols")
    (>&2 echo "rows: $rows")
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
		printf "%s\n" "${colname}" >> $tcount
		(>&2 echo "name: $colname")
      		(>&2 echo "--------------------------------------------------------------")
		cat $f >> $tcount
	    else
		touch $tfile1
		touch $tfile2
		colname=${f}
		#colname="${colname:2:${#colname}-1}"
		colname="${colname//$common_txt/}"
		colname="${colname%.*}"
		printf "%s\n" "${colname}" >> $tfile1
		(>&2 echo "name: $colname")
		(>&2 echo "--------------------------------------------------------------")
		cat $f | awk -v counts=$all_cols '{OFS="\t"; print $counts}' >> $tfile1
		paste $tcount $tfile1 >> $tfile2
		mv -f $tfile2 $tcount
		rm $tfile1
		#(>&2 echo $n)
	    fi
	    n=$(( $n + 1 ))
    fi 

done <<< "$DATA"

mv $tcount "${common_txt}gathered_counts.bed"
#mv $tcount test.bed



