INP_FILE=$1
OUTPUT_FOLDER=$2
COVERAGE=$3
SUPPORT_LARGE=$4
SUPPORT_IMP_LARGE=$5
SUPPORT_SMALL=$6
SUPPORT_SPLIT=$7

echo "Input file:" $INP_FILE
echo "Output folder:" $OUTPUT_FOLDER
echo "Coverage:" $COVERAGE
echo "Support Large:" $SUPPORT_LARGE
echo "Support Imprecise Large:" $SUPPORT_IMP_LARGE
echo "Support Small:" $SUPPORT_SMALL
echo "Support Split:" $SUPPORT_SPLIT

# Large
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp1_dels_large.txt" > $OUTPUT_FOLDER"tmp2_l0.txt"
cat $OUTPUT_FOLDER"tmp2_l0.txt" | awk -v SUP="$SUPPORT_LARGE" '{ if($5+$6+$7 >= SUP) { print }}' > $OUTPUT_FOLDER"tmp2_l1.txt"
samtools bedcov $OUTPUT_FOLDER"tmp2_l1.txt" $INP_FILE > $OUTPUT_FOLDER"tmp2_l2.txt"
awk -v COV="$COVERAGE" '{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= COV) {print $0, v[NR] }}' $OUTPUT_FOLDER"tmp2_l2.txt" > $OUTPUT_FOLDER"tmp2_l3.txt"
bedtools merge -i $OUTPUT_FOLDER"tmp2_l3.txt" > $OUTPUT_FOLDER"tmp2_large.bed"

# Imprecise Large
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp1_dels_large_imprecise.txt" > $OUTPUT_FOLDER"tmp2_il0.txt"
cat $OUTPUT_FOLDER"tmp2_il0.txt" | awk -v SUP="$SUPPORT_IMP_LARGE" '{ if($5+$6+$7 >= SUP) { print }}' > $OUTPUT_FOLDER"tmp2_il1.txt"
samtools bedcov $OUTPUT_FOLDER"tmp2_il1.txt" $INP_FILE > $OUTPUT_FOLDER"tmp2_il2.txt"
awk -v COV="$COVERAGE" '{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= COV) {print $0, v[NR] }}' $OUTPUT_FOLDER"tmp2_il2.txt" > $OUTPUT_FOLDER"tmp2_il3.txt"
bedtools merge -i $OUTPUT_FOLDER"tmp2_il3.txt" > $OUTPUT_FOLDER"tmp2_large_imprecise.bed"

# Split
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp1_dels_split.txt" > $OUTPUT_FOLDER"tmp2_sr0.txt"
cat $OUTPUT_FOLDER"tmp2_sr0.txt" | awk -v SUP="$SUPPORT_SPLIT" '{ if($5+$6+$7 >= SUP) { print }}' > $OUTPUT_FOLDER"tmp2_sr1.txt"
samtools bedcov $OUTPUT_FOLDER"tmp2_sr1.txt" $INP_FILE > $OUTPUT_FOLDER"tmp2_sr2.txt"
awk -v COV="$COVERAGE" '{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= COV) {print $0, v[NR] }}' $OUTPUT_FOLDER"tmp2_sr2.txt" > $OUTPUT_FOLDER"tmp2_sr3.txt"
bedtools merge -i $OUTPUT_FOLDER"tmp2_sr3.txt" > $OUTPUT_FOLDER"tmp2_split.bed"

# Small
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp1_dels_small.txt" > $OUTPUT_FOLDER"tmp2_s0.txt"
cat $OUTPUT_FOLDER"tmp2_s0.txt" | awk -v SUP="$SUPPORT_SMALL" '{ if($5+$6+$7 >= SUP) { print }}' > $OUTPUT_FOLDER"tmp2_s1.txt"
samtools bedcov $OUTPUT_FOLDER"tmp2_s1.txt" $INP_FILE > $OUTPUT_FOLDER"tmp2_s2.txt"
awk -v COV="$COVERAGE" '{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= COV) {print $0, v[NR] }}' $OUTPUT_FOLDER"tmp2_s2.txt" > $OUTPUT_FOLDER"tmp2_s3.txt"
bedtools merge -i $OUTPUT_FOLDER"tmp2_s3.txt" > $OUTPUT_FOLDER"tmp2_small.bed"

# Merge
cat $OUTPUT_FOLDER"tmp2_large.bed" $OUTPUT_FOLDER"tmp2_large_imprecise.bed" $OUTPUT_FOLDER"tmp2_small.bed" $OUTPUT_FOLDER"tmp2_split.bed" > $OUTPUT_FOLDER"tmp2_merge.bed"
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp2_merge.bed" > $OUTPUT_FOLDER"tmp2_merge_sort.bed"
bedtools merge -i $OUTPUT_FOLDER"tmp2_merge_sort.bed" > $OUTPUT_FOLDER"deletions.bed"

rm $OUTPUT_FOLDER"tmp2_l0.txt" $OUTPUT_FOLDER"tmp2_l1.txt" $OUTPUT_FOLDER"tmp2_l2.txt" $OUTPUT_FOLDER"tmp2_l3.txt" \
    $OUTPUT_FOLDER"tmp2_il0.txt" $OUTPUT_FOLDER"tmp2_il1.txt" $OUTPUT_FOLDER"tmp2_il2.txt" $OUTPUT_FOLDER"tmp2_il3.txt" \
    $OUTPUT_FOLDER"tmp2_s0.txt" $OUTPUT_FOLDER"tmp2_s1.txt" $OUTPUT_FOLDER"tmp2_s2.txt" $OUTPUT_FOLDER"tmp2_s3.txt" \
    $OUTPUT_FOLDER"tmp2_sr0.txt" $OUTPUT_FOLDER"tmp2_sr1.txt" $OUTPUT_FOLDER"tmp2_sr2.txt" $OUTPUT_FOLDER"tmp2_sr3.txt" \
    $OUTPUT_FOLDER"tmp2_merge.bed" $OUTPUT_FOLDER"tmp2_merge_sort.bed"