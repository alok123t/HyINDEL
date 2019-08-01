#!/bin/bash

# Absolute path to minia executable
# e.g., replace with "/path/to/minia"
if [[ -z $PATH_TO_MINIA ]]; then
    PATH_TO_MINIA="minia"
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
            esac
            shift
        done
    fi
}

function print_arguments() {
    echo "        Input file:" $INP_FILE >&2
    echo "        Output folder:" $OUTPUT_FOLDER >&2
}

parse_arguments $@

echo "[Step4] Insertions assembly start" >&2

print_arguments

$PATH_TO_SAMTOOLS view -@ 4 -f 12 -F 256 $INP_FILE | awk 'BEGIN{OFS=""}{($2==77)?f=1:f=2;printf("@%s/%d\n%s\n+\n%s\n", $1, f, $10, $11);}' >$OUTPUT_FOLDER"tmp/ins/orphans.fastq"
echo "        Orphans file:" $OUTPUT_FOLDER"tmp/ins/orphans.fastq" >&2

rm -f $OUTPUT_FOLDER"tmp/insertions.vcf"

for fName in "$OUTPUT_FOLDER"tmp/ins/*.fastq; do
    # Remove tmp files
    rm -f "${fName%.*}".h5 "${fName%.*}".unitigs* "${fName%.*}".contigs.fa
    # Call minia
    $PATH_TO_MINIA -verbose 0 -in "${fName}" -out "${fName%.*}" -kmer-size 31 >$OUTPUT_FOLDER"tmp/minia_stdout.txt" 2>$OUTPUT_FOLDER"tmp/minia_stderr.txt"
    # Remove tmp files
    rm "${fName%.*}".h5 "${fName%.*}".unitigs*
done

echo "[Step4] Insertions assembly end" >&2
