#!/bin/bash

function delete_tmp_files() {
    if [ "$1" = "all" ]; then
        rm -f $OUTPUT_FOLDER"tmp2_large.bed" $OUTPUT_FOLDER"tmp2_large_imprecise.bed" $OUTPUT_FOLDER"tmp2_small.bed"
    fi
    rm -f $OUTPUT_FOLDER"tmp2_l0.txt" $OUTPUT_FOLDER"tmp2_l1.txt" $OUTPUT_FOLDER"tmp2_l2.txt" $OUTPUT_FOLDER"tmp2_l3.txt" \
        $OUTPUT_FOLDER"tmp2_il0.txt" $OUTPUT_FOLDER"tmp2_il1.txt" $OUTPUT_FOLDER"tmp2_il2.txt" $OUTPUT_FOLDER"tmp2_il3.txt" \
        $OUTPUT_FOLDER"tmp2_s0.txt" $OUTPUT_FOLDER"tmp2_s1.txt" $OUTPUT_FOLDER"tmp2_s2.txt" $OUTPUT_FOLDER"tmp2_s3.txt" \
        $OUTPUT_FOLDER"tmp2_merge.bed" $OUTPUT_FOLDER"tmp2_merge_sort.bed"
}

ARGS=$(getopt i:o:c:l:m:s:q:h:: "$*")
eval set -- "$ARGS"

for arg; do
    case "$arg" in
    -i)
        INP_FILE=$2
        shift 2
        ;;
    -o)
        OUTPUT_FOLDER=$2
        shift 2
        ;;
    -c)
        COVERAGE=$2
        shift 2
        ;;
    -l)
        SUPPORT_LARGE=$2
        shift 2
        ;;
    -m)
        SUPPORT_IMP_LARGE=$2
        shift 2
        ;;
    -s)
        SUPPORT_SMALL=$2
        shift 2
        ;;
    -q)
        MAPPING_QUALITY=$2
        shift 2
        ;;
    -h)
        echo "Usage: bash /path/to/postProcess.sh -i /path/to/inputFile -o /path/to/outputFolder -c COV -l SUPPORT_LARGE -m SUPPORT_IMP_LARGE -s SUPPORT_SMALL -q MAPPING_QUALITY"
        exit 1
        ;;
    --)
        shift
        break
        ;;
    \?)
        echo "Invalid option: -$1"
        exit 0
        ;;
    esac
done

echo "Input file:" $INP_FILE
echo "Output folder:" $OUTPUT_FOLDER
echo "Coverage:" $COVERAGE
echo "Support Large:" $SUPPORT_LARGE
echo "Support Imprecise Large:" $SUPPORT_IMP_LARGE
echo "Support Small:" $SUPPORT_SMALL
echo "Mapping Quality:" $MAPPING_QUALITY

mkdir -p $OUTPUT_FOLDER

delete_tmp_files all

# Large
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp1_dels_large.txt" >$OUTPUT_FOLDER"tmp2_l0.txt"
cat $OUTPUT_FOLDER"tmp2_l0.txt" | awk -v SUP="$SUPPORT_LARGE" '{ if($5+$6+$7 >= SUP) { print }}' >$OUTPUT_FOLDER"tmp2_l1.txt"
samtools bedcov -Q $MAPPING_QUALITY $OUTPUT_FOLDER"tmp2_l1.txt" $INP_FILE >$OUTPUT_FOLDER"tmp2_l2.txt"
awk -v COV="$COVERAGE" '{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= COV) {print $0, v[NR] }}' $OUTPUT_FOLDER"tmp2_l2.txt" >$OUTPUT_FOLDER"tmp2_l3.txt"
bedtools merge -i $OUTPUT_FOLDER"tmp2_l3.txt" >$OUTPUT_FOLDER"tmp2_large.bed"

# Imprecise Large
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp1_dels_large_imprecise.txt" >$OUTPUT_FOLDER"tmp2_il0.txt"
cat $OUTPUT_FOLDER"tmp2_il0.txt" | awk -v SUP="$SUPPORT_IMP_LARGE" '{ if($5+$6+$7 >= SUP) { print }}' >$OUTPUT_FOLDER"tmp2_il1.txt"
samtools bedcov -Q $MAPPING_QUALITY $OUTPUT_FOLDER"tmp2_il1.txt" $INP_FILE >$OUTPUT_FOLDER"tmp2_il2.txt"
awk -v COV="$COVERAGE" '{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= COV) {print $0, v[NR] }}' $OUTPUT_FOLDER"tmp2_il2.txt" >$OUTPUT_FOLDER"tmp2_il3.txt"
bedtools merge -i $OUTPUT_FOLDER"tmp2_il3.txt" >$OUTPUT_FOLDER"tmp2_large_imprecise.bed"

# Small
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp1_dels_small.txt" >$OUTPUT_FOLDER"tmp2_s0.txt"
cat $OUTPUT_FOLDER"tmp2_s0.txt" | awk -v SUP="$SUPPORT_SMALL" '{ if($5+$6+$7 >= SUP) { print }}' >$OUTPUT_FOLDER"tmp2_s1.txt"
samtools bedcov -Q $MAPPING_QUALITY $OUTPUT_FOLDER"tmp2_s1.txt" $INP_FILE >$OUTPUT_FOLDER"tmp2_s2.txt"
awk -v COV="$COVERAGE" '{ a[NR]=$0; v[NR]=$NF/$4; if (v[NR] <= COV) {print $0, v[NR] }}' $OUTPUT_FOLDER"tmp2_s2.txt" >$OUTPUT_FOLDER"tmp2_s3.txt"
bedtools merge -i $OUTPUT_FOLDER"tmp2_s3.txt" >$OUTPUT_FOLDER"tmp2_small.bed"

# Merge
cat $OUTPUT_FOLDER"tmp2_large.bed" $OUTPUT_FOLDER"tmp2_large_imprecise.bed" $OUTPUT_FOLDER"tmp2_small.bed" >$OUTPUT_FOLDER"tmp2_merge.bed"
sort -k1,1 -k2,2n $OUTPUT_FOLDER"tmp2_merge.bed" >$OUTPUT_FOLDER"tmp2_merge_sort.bed"
bedtools merge -i $OUTPUT_FOLDER"tmp2_merge_sort.bed" >$OUTPUT_FOLDER"deletions.bed"

delete_tmp_files not_all
