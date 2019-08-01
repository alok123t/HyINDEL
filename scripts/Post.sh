#!/bin/bash

FLANK_SZ=1000
MIN_MAPPING_QUALITY=20

# Absolute path to samtools executable
# e.g., replace with "/path/to/samtools"
if [[ -z $PATH_TO_SAMTOOLS ]]; then
    PATH_TO_SAMTOOLS="samtools"
fi

function print_usage() {
    echo "Usage: $0 [-i /path/to/file] [-o /path/to/folder] [-c coverage]" >&2
    echo "        -h, --help  Help" >&2
    echo "        -i, --inp  Path to input file" >&2
    echo "        -o, --out  Path to output folder" >&2
    echo "        -c, --cov  Coverage" >&2
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
    echo "        Input file:" $INP_FILE >&2
    echo "        Output folder:" $OUT_FOLDER >&2
    echo "        Support Large:" $SUPPORT_LARGE >&2
    echo "        Support Imprecise Large:" $SUPPORT_IMP >&2
    echo "        Support Small:" $SUPPORT_SMALL >&2
    echo "        Coverage:" $COVERAGE >&2
}

function delete_tmp_files() {
    rm -rf $OUT_FOLDER"tmp/post/"
    rm -f $OUT_FOLDER"tmp/deletions.vcf"
}

function filter_support() {
    EVENTS=$OUT_FOLDER"tmp/dels/"$1
    FILTER_SUP_EVENTS=$OUT_FOLDER"tmp/post/0_"$1
    FILTER_DOUBLES=$OUT_FOLDER"tmp/post/0u_"$1
    cat $EVENTS | awk -v SUP="$2" '{ if($5+$6+$7 >= SUP) { print }}' >$FILTER_SUP_EVENTS
    sort -u $FILTER_SUP_EVENTS >$FILTER_DOUBLES
}

function filter_flanks_coverage() {
    IMP_FLAG=$2
    FILTER_DOUBLES=$OUT_FOLDER"tmp/post/0u_"$1
    SUP_FILE=$OUT_FOLDER"tmp/post/sup_"$1
    EVENT_FILE=$OUT_FOLDER"tmp/post/1_"$1
    A_FILE=$OUT_FOLDER"tmp/post/1a_"$1
    B_FILE=$OUT_FOLDER"tmp/post/1b_"$1
    COV_EVENT_FILE=$OUT_FOLDER"tmp/post/2_"$1
    COV_A_FILE=$OUT_FOLDER"tmp/post/2a_"$1
    COV_B_FILE=$OUT_FOLDER"tmp/post/2b_"$1
    COV_FINAL=$OUT_FOLDER"tmp/post/3_"$1
    FILTER_COV_EVENTS=$OUT_FOLDER"tmp/post/4_"$1
    OUT_FILE=$OUT_FOLDER"tmp/deletions.vcf"

    cat $FILTER_DOUBLES | awk '{ printf("%d\t%d\t%d\n", $5, $6, $7) }' >$SUP_FILE
    cat $FILTER_DOUBLES | awk -v FLANK="$FLANK_SZ" '{ printf("%s\t%d\t%d\n", $1, $2, $3) }' >$EVENT_FILE
    cat $FILTER_DOUBLES | awk -v FLANK="$FLANK_SZ" '{l=($2-FLANK);L=(l>0)?l:1; printf("%s\t%d\t%d\n", $1, L, $2) }' >$A_FILE
    cat $FILTER_DOUBLES | awk -v FLANK="$FLANK_SZ" '{ printf("%s\t%d\t%d\n", $1, $3, $3+FLANK) }' >$B_FILE

    $PATH_TO_SAMTOOLS bedcov -Q $MIN_MAPPING_QUALITY $EVENT_FILE $INP_FILE >$COV_EVENT_FILE
    $PATH_TO_SAMTOOLS bedcov -Q $MIN_MAPPING_QUALITY $A_FILE $INP_FILE >$COV_A_FILE
    $PATH_TO_SAMTOOLS bedcov -Q $MIN_MAPPING_QUALITY $B_FILE $INP_FILE >$COV_B_FILE

    paste $COV_EVENT_FILE $COV_A_FILE $COV_B_FILE $SUP_FILE >$COV_FINAL

    awk -v IMP_FLAG="$IMP_FLAG" '{ EV=$4/($3-$2+1); A=$8/1000; B=$12/1000; 
        printf("%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t", $1, $2, $3, $13, $14, $15, EV, A, B);
        if (A != 0 && B != 0) 
        {
            FL_A=EV/A; FL_B=EV/B; 
            if (FL_A >= 0.9 || FL_B >= 0.9) 
            {
                printf("0\t");
            } 
            else 
            {
                printf("1\t");
            }
        }
        else
        {
            printf("1\t");
        }
        printf("%d\n", IMP_FLAG);
    }' $COV_FINAL >$FILTER_COV_EVENTS

    awk '{ 
        printf("%s\t%d\t%d\t%s\t%s\t%s\t", $1, $2, "0", "N", "<DEL>", ".");
        if ($10==1)
        {
            printf("PASS\t");
        }
        else
        {
            printf(".\t");
        }
        printf("SVTYPE=DEL;SVLEN=%d;END=%d;SU=%d;PE=%d;SR=%d;SC=%d;COV=%.3f;COV_A=%.3f;COV_B=%.3f;", $3-$2+1, $3, $4+$5+$6, $4, $5, $6, $7, $8, $9);
        GT="./.";
        if ($11==1)
        {
            printf("IMPRECISE;");
        }
        else
        {
            if ($7 <= 0.2*$8 || $7 <= 0.2*$9)
            {
                GT="1/1";
            }
            else if ($7 <= 0.9*$8 && $7 <= 0.9*$9)
            {
                GT="0/1";
            }
            else
            {
                GT="./.";
            }
        }
        printf("\tGT:SU:PE:SR:SC\t%s:%d:%d:%d:%d\n", GT, $4+$5+$6, $4, $5, $6);
    }' $FILTER_COV_EVENTS >>$OUT_FILE
}

echo "[Step3] Postprocess start" >&2

parse_arguments $@

SUPPORT_LARGE=$(($COVERAGE / 3))
if [ "$SUPPORT_LARGE" -eq "0" ]; then
    SUPPORT_LARGE=1
fi
SUPPORT_SMALL=$(($COVERAGE / 6))
if [ "$SUPPORT_SMALL" -eq "0" ]; then
    SUPPORT_SMALL=1
fi
SUPPORT_IMP=$(($COVERAGE / 6))
if [ "$SUPPORT_IMP" -eq "0" ]; then
    SUPPORT_IMP=1
fi

print_arguments

delete_tmp_files

mkdir -p $OUT_FOLDER"tmp/post"

filter_support "small.txt" $SUPPORT_SMALL
filter_flanks_coverage "small.txt" 0

filter_support "large.txt" $SUPPORT_LARGE
filter_flanks_coverage "large.txt" 0

filter_support "large_imprecise.txt" $SUPPORT_IMP
filter_flanks_coverage "large_imprecise.txt" 1

echo "[Step3] Postprocess end" >&2
