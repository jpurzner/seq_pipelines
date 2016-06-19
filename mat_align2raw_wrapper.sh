
BAM_FILE=$1 
CUR_DIR="$(pwd)"
FILE_PREFIX="${BAM_FILE%.*}"
TAG_FILE="${CUR_DIR}/${FILE_PREFIX}.tagAlign"
UMAP_DIR="/tank/genomes/align2raw/umap/globalmap_k20tok54" 
OUT_FILE="${CUR_DIR}/${FILE_PREFIX}_norm.bedgraph"

MAT_COMMAND="align2rawsignal('-i=${TAG_FILE}', '-s=${mm9_chr_norandom}', '-u=${UMAP_DIR}', '-o=${OUT_FILE}', '-of=bg', '-mm=60', '-l=200');"

# check if tagAlign file already exists
if [ -a $TAG_FILE ];
then 
    echo "tagAlign file found"
else 
    echo "converting to tagAlign"
    samtools view -F 0x0204 ${BAM_FILE} | gawk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gawk 'BEGIN {OFS="\t"} {if ((($3-$2)>20) && (($3-$2)<54)) print $0; else if ((($3-$2)>54) && ($6=="+")) print $1, $2, $2+53, $4, $5, $6; else if ((($3-$2) > 54) && ($6=="-")) print $1, $3-53, $3, $4, $5, $6}'  | sort -k1,1 -k2,2n > ${TAG_FILE} 
echo "done" 
fi 

echo "starting align2rawsignal"
# create M file with function call and quit 
touch wiggler.m 
echo $MAT_COMMAND  >> wiggler.m 
echo "quit;" >> wiggler.m 

cat wiggler.m

matlab -nojvm -nodisplay -nosplash -nodesktop -r wiggler.m 
rm wiggler.m