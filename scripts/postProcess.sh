#!/bin/bash

FLANK_SZ=1000
MIN_MAPPING_QUALITY=20

# e.g., replace with "/path/to/samtools"
if [[ -z $PATH_TO_SAMTOOLS ]]; then
    PATH_TO_SAMTOOLS="samtools"
fi

function print_usage() {
    echo "Usage: $0 [-i /path/to/file] [-o /path/to/folder] [-c coverage]"
    echo "  -h, --help  Help"
    echo "  -i, --inp  Path to input file"
    echo "  -o, --out  Path to output folder"
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
    echo "Input file:" $INP_FILE
    echo "Output folder:" $OUT_FOLDER
    echo "Support Large:" $SUPPORT_LARGE
    echo "Support Imprecise Large:" $SUPPORT_IMP
    echo "Support Small:" $SUPPORT_SMALL
    echo "Coverage:" $COVERAGE
}

function delete_tmp_files() {
    rm -rf $OUT_FOLDER"tmp/post/"
}

function filter_support() {
    EVENTS=$OUT_FOLDER"tmp/dels/"$1
    FILTER_SUP_EVENTS=$OUT_FOLDER"tmp/post/0_"$1
    FILTER_DOUBLES=$OUT_FOLDER"tmp/post/0u_"$1
    cat $EVENTS | awk -v SUP="$2" '{ if($5+$6+$7 >= SUP) { print }}' >$FILTER_SUP_EVENTS
    sort -u $FILTER_SUP_EVENTS >$FILTER_DOUBLES
}

function filter_flanks_coverage() {
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
    OUT_FILE=$OUT_FOLDER"tmp/post/out.vcf"

    cat $FILTER_DOUBLES | awk '{ printf("%d\t%d\t%d\n", $5, $6, $7) }' >$SUP_FILE
    cat $FILTER_DOUBLES | awk -v FLANK="$FLANK_SZ" '{ printf("%s\t%d\t%d\n", $1, $2, $3) }' >$EVENT_FILE
    cat $FILTER_DOUBLES | awk -v FLANK="$FLANK_SZ" '{ printf("%s\t%d\t%d\n", $1, $2-FLANK, $2) }' >$A_FILE
    cat $FILTER_DOUBLES | awk -v FLANK="$FLANK_SZ" '{ printf("%s\t%d\t%d\n", $1, $3, $3+FLANK) }' >$B_FILE

    $PATH_TO_SAMTOOLS bedcov -Q $MIN_MAPPING_QUALITY $EVENT_FILE $INP_FILE >$COV_EVENT_FILE
    $PATH_TO_SAMTOOLS bedcov -Q $MIN_MAPPING_QUALITY $A_FILE $INP_FILE >$COV_A_FILE
    $PATH_TO_SAMTOOLS bedcov -Q $MIN_MAPPING_QUALITY $B_FILE $INP_FILE >$COV_B_FILE

    paste $COV_EVENT_FILE $COV_A_FILE $COV_B_FILE $SUP_FILE >$COV_FINAL

    awk '{ EV=$4/($3-$2+1); A=$8/1000; B=$12/1000; FL_A=EV/A; FL_B=EV/B; 
        printf("%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t", $1, $2, $3, $13, $14, $15, EV, A, B);
        if (FL_A >= 1.0 && FL_B >= 1.0) printf("0\n"); else printf("1\n"); 
        }' $COV_FINAL >$FILTER_COV_EVENTS

    cat $FILTER_COV_EVENTS | awk -v IMP_FLAG="$2" '{ id+=1;
        printf("%s\t%d\t%d\t%s\t%s\t%s\t", $1, $2, id, "N", "<DEL>", "."); 
        if ($10 == 1) printf("PASS\t"); else printf(".\t"); 
        printf("SVTYPE=DEL;SVLEN=%d;END=%d;SU=%d;PE=%d;SR=%d;SC=%d;COV=%.3f;COV_A=%.3f;COV_B=%.3f;", $3-$2+1, $3, $4+$5+$6, $4, $5, $6, $7, $8, $9);
        if ($IMP_FLAG == 1) printf("IMPRECISE;");
        printf("\tSU:PE:SR:SC\t%d:%d:%d:%d\n", $4+$5+$6, $4, $5, $6);
        }' >>$OUT_FILE
}

echo "[Step5] Postprocess start"

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

OUT_FILE=$OUT_FOLDER"tmp/post/out.vcf"
OUT_SORT_FILE=$OUT_FOLDER"tmp/post/out_sort.vcf"
FINAL_FILE=$OUT_FOLDER"output.vcf"

sort -k1,1 -k2,2n $OUT_FILE >$OUT_SORT_FILE
cat $OUT_FOLDER"tmp/pre/header.vcf" >$FINAL_FILE
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample" >>$FINAL_FILE
cat $OUT_SORT_FILE >>$FINAL_FILE

echo "[Step5] Postprocess end"
