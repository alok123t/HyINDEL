#!/bin/bash

# Absolute path to minia3 executable
# e.g., replace with "/path/to/minia"
if [[ -z $PATH_TO_MINIA ]]; then
    PATH_TO_MINIA="minia"
fi

# Absolute path to minimap2 executable
# e.g., replace with "/path/to/minia"
if [[ -z $PATH_TO_MINIMAP ]]; then
    PATH_TO_MINIMAP="minimap2"
fi

# Absolute path to samtools executable
# e.g., replace with "/path/to/samtools"
if [[ -z $PATH_TO_SAMTOOLS ]]; then
    PATH_TO_SAMTOOLS="samtools"
fi

function print_usage() {
    echo "Usage: $0 [-i /path/to/file] [-o /path/to/folder]" >&2
    echo "        -h, --help  Help" >&2
    echo "        -i, --inp  Path to input file" >&2
    echo "        -o, --out  Path to output folder" >&2
}

function parse_arguments() {
    if [[ -z $1 ]]; then
        print_usage
        exit
    else
        while [ "$1" != "" ]; do
            case $1 in
            -h | --help)
                print_usage
                exit
                ;;
            -i | --inp)
                shift
                INP_FILE=$1
                ;;
            -o | --out)
                shift
                OUTPUT_FOLDER=$1
                # Add / at end if not present
                LAST_CHAR_OUT_FOLDER="${OUTPUT_FOLDER: -1}"
                if [ $LAST_CHAR_OUT_FOLDER != "/" ]; then
                    OUTPUT_FOLDER="${OUTPUT_FOLDER}/"
                fi
                ;;
            -r | --ref)
                shift
                REF_FILE=$1
                ;;
            esac
            shift
        done
    fi
}

function print_arguments() {
    echo "        Input file:" $INP_FILE >&2
    echo "        Output folder:" $OUTPUT_FOLDER >&2
    echo "        Reference index:" $REF_FILE >&2
}

parse_arguments $@

echo "[Step4] Assembly start" >&2

print_arguments

rm -f $OUTPUT_FOLDER"tmp/insertions.vcf"

# Append orphan reads
$PATH_TO_SAMTOOLS view -@ 4 -f 12 -F 256 $INP_FILE | awk 'BEGIN{OFS=""}{($2==77)?f=1:f=2;printf("@%s/%d\n%s\n+\n%s\n", $1, f, $10, $11);}' >>$OUTPUT_FOLDER"tmp/ins/reads.fastq"

# Assemble reads using minia
$PATH_TO_MINIA -verbose 0 -in $OUTPUT_FOLDER"tmp/ins/reads.fastq" -out $OUTPUT_FOLDER"tmp/ins/31" -kmer-size 31 >$OUTPUT_FOLDER"tmp/minia_stdout.txt" 2>$OUTPUT_FOLDER"tmp/minia_stderr.txt"
rm $OUTPUT_FOLDER"tmp/ins/reads.fastq"
rm $OUTPUT_FOLDER"tmp/ins/31".unitigs* $OUTPUT_FOLDER"tmp/ins/31".h5

# Align contigs to reference using minimap
$PATH_TO_MINIMAP -a $REF_FILE $OUTPUT_FOLDER"tmp/ins/31".contigs.fa >$OUTPUT_FOLDER"tmp/ins/31_contigs.sam"

# Sort and convert sam to bam
$PATH_TO_SAMTOOLS view -bS $OUTPUT_FOLDER"tmp/ins/31_contigs.sam" >$OUTPUT_FOLDER"tmp/ins/31_contigs.bam"
$PATH_TO_SAMTOOLS sort $OUTPUT_FOLDER"tmp/ins/31_contigs.bam" -o $OUTPUT_FOLDER"tmp/ins/31_contigs_sort.bam"
# Index bam file
$PATH_TO_SAMTOOLS index $OUTPUT_FOLDER"tmp/ins/31_contigs_sort.bam"
rm -f $OUTPUT_FOLDER"tmp/ins/31_contigs.sam" $OUTPUT_FOLDER"tmp/ins/31_contigs.bam"

echo "[Step4] Assembly end" >&2
