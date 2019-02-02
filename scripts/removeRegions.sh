WINDOW_SIZE=1000

INP_FILE=$1
OUTPUT_FOLDER=$2

echo "Input file:" $INP_FILE
echo "Output folder:" $OUTPUT_FOLDER

# # Get chromosome sizes from header
# samtools idxstats $INP_FILE | cut -f 1,2 > $OUTPUT_FOLDER"tmp0_chr_sizes.bed"

# # Divide each chromosome into windows/bins
# bedtools makewindows -g $OUTPUT_FOLDER"tmp0_chr_sizes.bed" -w $WINDOW_SIZE > $OUTPUT_FOLDER"tmp1_chr_windows.bed"

# # Number of reads in each window
# samtools bedcov $OUTPUT_FOLDER"tmp1_chr_windows.bed" $INP_FILE > $OUTPUT_FOLDER"tmp2_chr_windows_count.bed"

# # Get coverage for each window
# awk -f div-cov.awk $OUTPUT_FOLDER"tmp2_chr_windows_count.bed" > $OUTPUT_FOLDER"tmp3_chr_windows_cov.bed"

# Reverse sort file by coverage 
# sort -k 4 -n $OUTPUT_FOLDER"tmp0_cov.regions.bed" > $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed"

# Calculate median of coverage

# MEDIAN=$(awk '{a[i++]=$4} END {x=int((i+1)/2); if (x<(i+1)/2) print(int(a[x-1]+a[x])/2); else print int(a[x-1]);}' $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed")
# echo "Median:" $MEDIAN
# EXCLUDE_CUTOFF=$(($MEDIAN*3))
# echo "Exclude Cutoff:" $EXCLUDE_CUTOFF

EXCLUDE_CUTOFF=99

cat $OUTPUT_FOLDER"tmp0_cov.regions.sort.bed" | awk -v EXC="$EXCLUDE_CUTOFF" '{if ($4 >= EXC) {printf("%d\t%d\t%d\n", $1,$2,$3)}}' > $OUTPUT_FOLDER"remove_chr_windows.bed"

# # Number of rows to be removed
# REMOVE_LINES=$(wc -l $OUTPUT_FOLDER"tmp4_chr_windows_cov_rsort.bed" | awk '{printf "%d", $1*0.005}')

# # Final output file containing windows to be excluded
# head -n $REMOVE_LINES $OUTPUT_FOLDER"tmp4_chr_windows_cov_rsort.bed" > $OUTPUT_FOLDER"remove_chr_windows.bed"

# # Remove tmp files
# rm $OUTPUT_FOLDER"tmp0_chr_sizes.bed" $OUTPUT_FOLDER"tmp1_chr_windows.bed" $OUTPUT_FOLDER"tmp2_chr_windows_count.bed" $OUTPUT_FOLDER"tmp3_chr_windows_cov.bed"
