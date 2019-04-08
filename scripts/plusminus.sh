#!/bin/bash

# Absolute path to executable e.g., /usr/local/bin/mosdepth
export PATH_TO_MOSDEPTH=""

export PATH_TO_MINIA=""

function print_usage() {
    echo "Usage: $0 [-i /path/to/file] [-o /path/to/folder] [-s insert size] [-d standard deviation] [-l read length] [-c coverage]"
    echo "  -h, --help  Help"
    echo "  -i, --inp  Path to input file"
    echo "  -o, --out  Path to output folder"
    echo "  -s, --insSz  Median Insert size"
    echo "  -d, --stdDev  Median Standard deviation"
    echo "  -l, --readLen  Read length"
    echo "  -c, --cov  Coverage"
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
            esac
            shift
        done
    fi
}

function print_arguments() {
    INVALID_ARG=0
    if [[ -z $INP_FILE ]]; then
        echo "ERROR: -i missing"
        INVALID_ARG=1
    else
        if [[ -e $INP_FILE ]]; then
            echo "Input file:" $INP_FILE
        else
            echo $INP_FILE "does not exist"
            exit
        fi
    fi
    if [[ -z $OUT_FOLDER ]]; then
        echo "ERROR: -o missing"
        INVALID_ARG=1
    else
        if [[ -d $OUT_FOLDER ]]; then
            echo "Output folder:" $OUT_FOLDER
        else
            # Create folder if it doesn't exist
            mkdir -p $OUT_FOLDER
            echo $OUT_FOLDER "created"
        fi
    fi
    if [[ -z $INS_SZ ]]; then
        echo "ERROR: -s missing"
        INVALID_ARG=1
    else
        if [[ $INS_SZ -gt 0 ]]; then
            echo "Insert size:" $INS_SZ
        else
            echo $INS_SZ "is invalid"
            exit
        fi
    fi
    if [[ -z $STD_DEV ]]; then
        echo "ERROR: -d missing"
        INVALID_ARG=1
    else
        if [[ $STD_DEV -gt 0 ]]; then
            echo "Standard deviation:" $STD_DEV
        else
            echo $STD_DEV "is invalid"
            exit
        fi
    fi
    if [[ -z $INS_SZ ]]; then
        echo "ERROR: -l missing"
        INVALID_ARG=1
    else
        if [[ $READ_LEN -gt 0 ]]; then
            echo "Read length:" $READ_LEN
        else
            echo $READ_LEN "is invalid"
            exit
        fi
    fi
    if [[ -z $COVERAGE ]]; then
        echo "ERROR: -c missing"
        INVALID_ARG=1
    else
        if [[ $COVERAGE -gt 0 ]]; then
            echo "Coverage:" $COVERAGE
        else
            echo $COVERAGE "is invalid"
            exit
        fi
    fi

    if [[ $INVALID_ARG -eq 1 ]]; then
        print_usage
        exit
    fi
}

parse_arguments $@

print_arguments

# Change to directory containing this script
cd $(dirname $0)

./preProcess -i $INP_FILE -o $OUT_FOLDER

./indel -i $INP_FILE -o $OUT_FOLDER -s $INS_SZ -d $STD_DEV -l $READ_LEN

./assembleInsertions -o $OUT_FOLDER

./ins -o $OUT_FOLDER

./postProcess -i $INP_FILE -o $OUT_FOLDER -c $COVERAGE -l 25 -m 15 -s 15 -q 20
