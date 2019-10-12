#!/bin/bash

# Absolute path to mosdepth executable
# e.g., replace with "/path/to/mosdepth"
if [[ -z $PATH_TO_MOSDEPTH ]]; then
    PATH_TO_MOSDEPTH="mosdepth"
fi

WINDOW_SIZE=10000000
OVERLAP_SIZE=20000
SLIDING_SIZE=$((WINDOW_SIZE - OVERLAP_SIZE))

REMOVE_WINDOW_SIZE=1000
MOS_THREADS=4

SKIP_CALC=0

function print_usage() {
    echo "Usage: $0 [-i /path/to/file] [-o /path/to/folder]" >&2
    echo "  -h, --help  Help" >&2
    echo "  -i, --inp  Path to input file" >&2
    echo "  -o, --out  Path to output folder" >&2
}

function print_vcf_header() {
    VCF=$OUT_FOLDER"tmp/pre/header.vcf"
    echo "##fileformat=VCFv4.2" >$VCF
    echo "##source=HyINDEL" >>$VCF
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">" >>$VCF
    echo "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">" >>$VCF
    echo "##INFO=<ID=SEQ,Number=.,Type=String,Description="Novel sequence insertion">" >>$VCF
    echo "##INFO=<ID=SEQUP,Number=.,Type=String,Description="Partial upstream sequence insertion">" >>$VCF
    echo "##INFO=<ID=SEQDOWN,Number=.,Type=String,Description="Partial downstream sequence insertion">" >>$VCF
    echo "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">" >>$VCF
    echo "##INFO=<ID=SUP,Number=.,Type=Integer,Description="Total number of reads supporting this variant">" >>$VCF
    echo "##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting this variant">" >>$VCF
    echo "##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting this variant">" >>$VCF
    echo "##INFO=<ID=SC,Number=.,Type=Integer,Description="Number of softclip reads supporting this variant">" >>$VCF
    echo "##INFO=<ID=COV,Number=.,Type=Float,Description="Coverage in between breakpoints">" >>$VCF
    echo "##ALT=<ID=DEL,Description="Deletion">" >>$VCF
    echo "##ALT=<ID=INS,Description="Insertion">" >>$VCF

    TMP_HEADER=$OUT_FOLDER"tmp/pre/chr_sizes.bed"
    cat $TMP_HEADER | awk '{ cid=$2;clen=$3;split(cid, cidAr, ":");split(clen, clenAr, ":");if(cidAr[1]=="SN") print("##contig=<ID=" cidAr[2] ",length=" clenAr[2]) ">"}' >>$VCF
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
            --skip)
                shift
                SKIP_CALC=$1
                ;;
            esac
            shift
        done
    fi
}

function print_arguments() {
    echo "        Input file:" $INP_FILE >&2
    echo "        Output folder:" $OUT_FOLDER >&2
}

echo "[Step1] Preprocess start" >&2

parse_arguments $@

print_arguments

# Get chromosome sizes from header
$(dirname "$0")/bamtools header -in $INP_FILE >$OUT_FOLDER"tmp/pre/chr_sizes.bed"

# Print vcf header
print_vcf_header

# Divide each chromosome into overlapping windows/bins for variant processing
cat $OUT_FOLDER"tmp/pre/chr_sizes.bed" | awk -v WIN_SZ="$WINDOW_SIZE" -v SLD_SZ="$SLIDING_SIZE" '{if ($1 == "@SQ") {split($2,a,":");split($3,b,":"); for (i = 0; i <= b[2]; i+=SLD_SZ) {en=(i+WIN_SZ>b[2])?b[2]:i+WIN_SZ;printf("%s\t%d\t%d\n",a[2],i,en);}}}' >$OUT_FOLDER"tmp/pre/chr_windows.bed"

if [ "$SKIP_CALC" -eq "1" ]; then
    touch $OUT_FOLDER"tmp/pre/remove_chr_windows.bed"
    echo "[Step1] Preprocess skip coverage calculation" >&2
    exit
fi

# Calculate coverage for each removal window
$PATH_TO_MOSDEPTH -n --fast-mode -t $MOS_THREADS --by $REMOVE_WINDOW_SIZE $OUT_FOLDER"tmp/pre/cov" $INP_FILE
gunzip -c $OUT_FOLDER"tmp/pre/cov.regions.bed.gz" >$OUT_FOLDER"tmp/pre/cov.regions.bed"

# Sort file by coverage
sort -k 4 -n $OUT_FOLDER"tmp/pre/cov.regions.bed" >$OUT_FOLDER"tmp/pre/cov.regions.sort.bed"

# Calculate median of coverage
MEDIAN=$(awk '{if ($4 != 0) a[i++]=$4} END {x=int((i+1)/2); if (x<(i+1)/2) print(int((a[x-1]+a[x])/2)); else print int(a[x-1]);}' $OUT_FOLDER"tmp/pre/cov.regions.sort.bed")
echo "        Median:" $MEDIAN >&2
EXCLUDE_CUTOFF=$(($MEDIAN * 3))
echo "        Exclude Cutoff:" $EXCLUDE_CUTOFF >&2

# Filter windows by exclude cutoff
cat $OUT_FOLDER"tmp/pre/cov.regions.sort.bed" | awk -v EXC="$EXCLUDE_CUTOFF" '{if ($4 >= EXC) {printf("%s\t%d\t%d\n", $1,$2,$3)}}' >$OUT_FOLDER"tmp/pre/remove_chr_windows.bed"

# Delete large tmp files
rm $OUT_FOLDER"tmp/pre/cov.regions.bed" $OUT_FOLDER"tmp/pre/cov.regions.sort.bed" $OUT_FOLDER"tmp/pre/cov.regions.bed.gz" $OUT_FOLDER"tmp/pre/cov.mosdepth.region.dist.txt" $OUT_FOLDER"tmp/pre/cov.mosdepth.global.dist.txt"

echo "[Step1] Preprocess end" >&2
