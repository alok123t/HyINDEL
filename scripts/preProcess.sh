#!/bin/bash

WINDOW_SIZE=10000000
OVERLAP_SIZE=20000
SLIDING_SIZE=$((WINDOW_SIZE - OVERLAP_SIZE))

REMOVE_WINDOW_SIZE=1000
MOS_THREADS=4

ARGS=$(getopt i:o:h::m: "$*")
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
    -m)
        USE_MOSDEPTH=$2
        shift 2
        ;;
    -h)
        echo "Usage: bash /path/to/preProcess.sh -i /path/to/inputFile -o /path/to/outputFolder -m 1"
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

LAST_CHAR_OUT_FOLDER="${OUTPUT_FOLDER: -1}"
if [ $LAST_CHAR_OUT_FOLDER != "/" ]; then
    OUTPUT_FOLDER="${OUTPUT_FOLDER}/"
fi

echo "Input file:" $INP_FILE
echo "Output folder:" $OUTPUT_FOLDER

mkdir -p $OUTPUT_FOLDER/tmp_ins

# Extract unmapped reads
bamtools filter -in $INP_FILE -out $OUTPUT_FOLDER"unmapped.bam" -script "$(dirname "$0")"/filter_unmapped.json

# Convert unmapped bam to fastq
bamtools convert -in $OUTPUT_FOLDER"unmapped.bam" -out $OUTPUT_FOLDER"unmapped.fastq" -format fastq
rm $OUTPUT_FOLDER"unmapped.bam"
echo "Created" $OUTPUT_FOLDER"unmapped.fastq"

# Get chromosome sizes from header
bamtools header -in $INP_FILE >$OUTPUT_FOLDER"tmp0_chr_sizes.bed"

# Divide each chromosome into overlapping windows/bins for variant processing
cat $OUTPUT_FOLDER"tmp0_chr_sizes.bed" | awk -v WIN_SZ="$WINDOW_SIZE" -v SLD_SZ="$SLIDING_SIZE" '{if ($1 == "@SQ") {split($2,a,":");split($3,b,":"); for (i = 0; i <= b[2]; i+=SLD_SZ) {en=(i+WIN_SZ>b[2])?b[2]:i+WIN_SZ;printf("%s\t%d\t%d\n",a[2],i,en);}}}' >$OUTPUT_FOLDER"chr_windows.bed"

# Calculate coverage for each removal window
if [ $USE_MOSDEPTH = 1 ]; then
    echo "Using mosdepth for coverage calculation"
    mosdepth -n --fast-mode -t $MOS_THREADS --by $REMOVE_WINDOW_SIZE $OUTPUT_FOLDER"cov" $INP_FILE
    gunzip $OUTPUT_FOLDER"cov.regions.bed.gz"
else
    echo "Using samtools for coverage calculation"
    cat $OUTPUT_FOLDER"tmp0_chr_sizes.bed" | awk -v REM_SZ="$REMOVE_WINDOW_SIZE" '{if ($1 == "@SQ") {split($2,a,":");split($3,b,":"); for (i = 0; i <= b[2]; i+=REM_SZ) {en=(i+REM_SZ>b[2])?b[2]:i+REM_SZ;printf("%s\t%d\t%d\n",a[2],i,en);}}}' >$OUTPUT_FOLDER"tmp0_chr_rem_windows.bed"
    samtools bedcov $OUTPUT_FOLDER"tmp0_chr_rem_windows.bed" $INP_FILE >$OUTPUT_FOLDER"tmp0_cov.regions.bed"
    awk -v SZ="$REMOVE_WINDOW_SIZE" '{ cov=$4/SZ; {print $1,$2,$3,cov}}' $OUTPUT_FOLDER"tmp0_cov.regions.bed" >$OUTPUT_FOLDER"cov.regions.bed"
fi

# Sort file by coverage
sort -k 4 -n $OUTPUT_FOLDER"cov.regions.bed" >$OUTPUT_FOLDER"tmp0_cov.regions.sort.bed"

# Calculate median of coverage
MEDIAN=$(awk '{if ($4 != 0) a[i++]=$4} END {x=int((i+1)/2); if (x<(i+1)/2) print(int((a[x-1]+a[x])/2)); else print int(a[x-1]);}' $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed")
echo "Median:" $MEDIAN
EXCLUDE_CUTOFF=$(($MEDIAN * 2))
echo "Exclude Cutoff:" $EXCLUDE_CUTOFF

# Filter windows by exclude cutoff
cat $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed" | awk -v EXC="$EXCLUDE_CUTOFF" '{if ($4 >= EXC) {printf("%s\t%d\t%d\n", $1,$2,$3)}}' >$OUTPUT_FOLDER"remove_chr_windows.bed"

# Remove tmp files
rm -f $OUTPUT_FOLDER"tmp0_chr_sizes.bed" $OUTPUT_FOLDER"tmp0_chr_rem_windows.bed" $OUTPUT_FOLDER"tmp0_cov.regions.bed" $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed" $OUTPUT_FOLDER"cov.regions.bed"
