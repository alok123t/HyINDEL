INP_FILE=$1
OUTPUT_FOLDER=$2

WINDOW_SIZE=1000000
SLIDING_SIZE=800000
REMOVE_WINDOW_SIZE=1000
MOS_THREADS=4

echo "Input file:" $INP_FILE
echo "Output folder:" $OUTPUT_FOLDER

# Get chromosome sizes from header
samtools idxstats $INP_FILE | cut -f 1,2 > $OUTPUT_FOLDER"tmp0_chr_sizes.bed"

# Divide each chromosome into overlapping windows/bins for variant processing
bedtools makewindows -g $OUTPUT_FOLDER"tmp0_chr_sizes.bed" -w $WINDOW_SIZE -s $SLIDING_SIZE > $OUTPUT_FOLDER"chr_windows.bed"

# Calculate coverage for each removal window
mosdepth -n --fast-mode -t $MOS_THREADS --by $REMOVE_WINDOW_SIZE $OUTPUT_FOLDER"cov" $INP_FILE

# Extract gzip file
gunzip $OUTPUT_FOLDER"cov.regions.bed.gz"

# Sort file by coverage
sort -k 4 -n $OUTPUT_FOLDER"cov.regions.bed" > $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed"

# Calculate median of coverage
MEDIAN=$(awk '{a[i++]=$4} END {x=int((i+1)/2); if (x<(i+1)/2) print(int(a[x-1]+a[x])/2); else print int(a[x-1]);}' $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed")
echo "Median:" $MEDIAN
EXCLUDE_CUTOFF=$(($MEDIAN*2))
echo "Exclude Cutoff:" $EXCLUDE_CUTOFF

# Filter windows by exclude cutoff
cat $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed" | awk -v EXC="$EXCLUDE_CUTOFF" '{if ($4 >= EXC) {printf("%s\t%d\t%d\n", $1,$2,$3)}}' > $OUTPUT_FOLDER"remove_chr_windows.bed"

# Remove tmp files
rm $OUTPUT_FOLDER"tmp0_chr_sizes.bed" $OUTPUT_FOLDER"tmp0_remove_chr_windows.bed" $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed"