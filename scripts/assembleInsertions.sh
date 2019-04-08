#!/bin/bash

# Absolute path to MINIA executable
if [[ -z $PATH_TO_MINIA ]]; then
    PATH_TO_MINIA="/Users/alok/Tools/minia/build/bin/minia"
fi

ARGS=$(getopt o:h:: "$*")
eval set -- "$ARGS"

for arg; do
    case "$arg" in
    -o)
        OUTPUT_FOLDER=$2
        shift 2
        ;;
    -h)
        echo "Usage: bash /path/to/assembleInsertions.sh -o /path/to/outputFolder"
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

echo "[Step3] Assemble insertions start"
echo "Output folder:" $OUTPUT_FOLDER

for fName in "$OUTPUT_FOLDER"tmp/ins/*.fastq; do
    $PATH_TO_MINIA -verbose 0 -in "${fName}" -out "${fName%.*}" >>$OUTPUT_FOLDER"tmp/minia_stdout.txt" 2>>$OUTPUT_FOLDER"tmp/minia_stderr.txt"
done

echo "[Step3] Assemble insertions end"
