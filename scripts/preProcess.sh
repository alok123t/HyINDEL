INP_FILE=$1
OUTPUT_FOLDER=$2

WINDOW_SIZE=1000000

echo "Input file:" $INP_FILE
echo "Output folder:" $OUTPUT_FOLDER

samtools idxstats $INP_FILE | cut -f 1,2 > $OUTPUT_FOLDER"chr_sizes.bed"
bedtools makewindows -g $OUTPUT_FOLDER"chr_sizes.bed" -w $WINDOW_SIZE > $OUTPUT_FOLDER"chr_windows.bed"
