#!/bin/bash

THREADS=1
SKIP_PRE=0

function print_usage() {
    echo "Usage: $0 [-i /path/to/file] [-o /path/to/folder] [-r /path/to/file] [-s insert size] [-d standard deviation] [-l read length] [-c coverage] [-t threads]" >&2
    echo "  -h, --help  Help" >&2
    echo "  -i, --inp  Path to input file" >&2
    echo "  -o, --out  Path to output folder" >&2
    echo "  -r, --ref  Path to reference file" >&2
    echo "  -s, --insSz  Median Insert size" >&2
    echo "  -d, --stdDev  Median Standard deviation" >&2
    echo "  -l, --readLen  Read length" >&2
    echo "  -c, --cov  Coverage" >&2
    echo "  -t, --threads  Threads" >&2
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
                OUT_FOLDER=$1
                # Add / at end if not present
                LAST_CHAR_OUT_FOLDER="${OUT_FOLDER: -1}"
                if [ $LAST_CHAR_OUT_FOLDER != "/" ]; then
                    OUT_FOLDER="${OUT_FOLDER}/"
                fi
                ;;
            -r | --ref)
                shift
                REF_FILE=$1
                ;;
            -s | --insSz)
                shift
                INS_SZ=$1
                ;;
            -d | --stdDev)
                shift
                STD_DEV=$1
                ;;
            -l | --readLen)
                shift
                READ_LEN=$1
                ;;
            -c | --cov)
                shift
                COVERAGE=$1
                ;;
            -t | --threads)
                shift
                THREADS=$1
                ;;
            --skip)
                shift
                SKIP_PRE=$1
                ;;
            esac
            shift
        done
    fi
}

function print_arguments() {
    INVALID_ARG=0
    if [[ -z $INP_FILE ]]; then
        echo "ERROR: -i missing" >&2
        INVALID_ARG=1
    else
        if [[ -e $INP_FILE ]]; then
            echo "        Input file:" $INP_FILE >&2
        else
            echo "ERROR: " $INP_FILE "does not exist" >&2
            exit
        fi
    fi
    if [[ -z $OUT_FOLDER ]]; then
        echo "ERROR: -o missing" >&2
        INVALID_ARG=1
    else
        if [[ -d $OUT_FOLDER ]]; then
            echo "        Output folder:" $OUT_FOLDER >&2
        else
            # Create folder if it doesn't exist
            mkdir -p $OUT_FOLDER
            echo "        Output folder:" $OUT_FOLDER "created" >&2
        fi
    fi
    if [[ -z $REF_FILE ]]; then
        echo "ERROR: -r missing" >&2
        INVALID_ARG=1
    else
        if [[ -e $REF_FILE ]]; then
            echo "        Reference file:" $REF_FILE >&2
        else
            echo "ERROR: " $REF_FILE "does not exist" >&2
            exit
        fi
    fi
    if [[ -z $INS_SZ ]]; then
        echo "ERROR: -s missing" >&2
        INVALID_ARG=1
    else
        if [[ $INS_SZ -gt 0 ]]; then
            echo "        Insert size:" $INS_SZ >&2
        else
            echo "ERROR: Insert size" $INS_SZ "is invalid" >&2
            exit
        fi
    fi
    if [[ -z $STD_DEV ]]; then
        echo "ERROR: -d missing" >&2
        INVALID_ARG=1
    else
        if [[ $STD_DEV -gt 0 ]]; then
            echo "        Standard deviation:" $STD_DEV >&2
        else
            echo "ERROR: Standard deviation" $STD_DEV "is invalid" >&2
            exit
        fi
    fi
    if [[ -z $INS_SZ ]]; then
        echo "ERROR: -l missing" >&2
        INVALID_ARG=1
    else
        if [[ $READ_LEN -gt 0 ]]; then
            echo "        Read length:" $READ_LEN >&2
        else
            echo "ERROR: Read length" $READ_LEN "is invalid" >&2
            exit
        fi
    fi
    if [[ -z $COVERAGE ]]; then
        echo "ERROR: -c missing" >&2
        INVALID_ARG=1
    else
        if [[ $COVERAGE -gt 0 ]]; then
            echo "        Coverage:" $COVERAGE >&2
        else
            echo "ERROR: Coverage" $COVERAGE "is invalid" >&2
            exit
        fi
    fi
    if [[ $THREADS -gt 0 ]]; then
        echo "        Threads:" $THREADS >&2
    else
        echo "ERROR: Threads" $THREADS "is invalid" >&2
        exit
    fi

    if [[ $INVALID_ARG -eq 1 ]]; then
        print_usage
        exit
    fi
}

parse_arguments $@

print_arguments

mkdir -p $OUT_FOLDER/tmp/pre $OUT_FOLDER/tmp/dels $OUT_FOLDER/tmp/ins $OUT_FOLDER/tmp/post

$(dirname $0)/HyINDEL_pre -i $INP_FILE -o $OUT_FOLDER --skip $SKIP_PRE

$(dirname $0)/HyINDEL_dels -i $INP_FILE -o $OUT_FOLDER -s $INS_SZ -d $STD_DEV -l $READ_LEN -c $COVERAGE -t $THREADS

$(dirname $0)/HyINDEL_post -i $INP_FILE -o $OUT_FOLDER -c $COVERAGE

$(dirname $0)/HyINDEL_assembly -i $INP_FILE -o $OUT_FOLDER -r $REF_FILE

$(dirname $0)/HyINDEL_ins -i $OUT_FOLDER"tmp/ins/31_contigs_sort.bam" -o $OUT_FOLDER

DELS_FILE=$OUT_FOLDER"tmp/deletions.vcf"
INS_FILE=$OUT_FOLDER"tmp/insertions.vcf"
OUT_FILE=$OUT_FOLDER"tmp/out_1.vcf"
OUT_SORT_FILE=$OUT_FOLDER"tmp/out_2.vcf"
OUT_ID_FILE=$OUT_FOLDER"tmp/out_3.vcf"
FINAL_FILE=$OUT_FOLDER"output.vcf"

cat $DELS_FILE $INS_FILE >$OUT_FILE

sort -k1,1 -k2,2n $OUT_FILE >$OUT_SORT_FILE

cat $OUT_SORT_FILE | awk '{ co+=1; printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, co, $4, $5, $6, $7, $8, $9, $10);}' >$OUT_ID_FILE

cat $OUT_FOLDER"tmp/pre/header.vcf" >$FINAL_FILE
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample" >>$FINAL_FILE
cat $OUT_ID_FILE >>$FINAL_FILE

rm $OUT_FILE $OUT_SORT_FILE $OUT_ID_FILE
