#!/bin/bash
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -N 1
#SBATCH -n 22
#SBATCH -p short
#SBATCH -w node48
#SBATCH --time=0-06:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-type=END

Time="/usr/bin/time"
Bwa="/home/iiit/alok.t/miniconda2/bin/bwa"
Samtools="/home/iiit/alok.t/miniconda2/bin/samtools"
Fastqc="/home/iiit/alok.t/Tools/FastQC/fastqc"
Bamqc="/home/iiit/alok.t/Tools/BamQC/bin/bamqc"
Picard="/home/iiit/alok.t/Tools/picard.jar"
Tiddit="/home/iiit/alok.t/Tools/TIDDIT-TIDDIT-2.6.0/TIDDIT.py"
SoftSV="/home/iiit/alok.t/Tools/SoftSV_1.4.2/SoftSV"
Pamir="/home/iiit/alok.t/Tools/pamir-master/pamir.py"
HyINDEL="/home/iiit/alok.t/install/HyINDEL/bin/HyINDEL"
GROM="/home/iiit/alok.t/Tools/GROM/dist/GROM"
Popins="/home/iiit/alok.t/Tools/popins/popins"

DIR="/scratch/alok/Platinum/"
mkdir -p $DIR
mkdir -p $DIR"Reference/"
mkdir -p $DIR"Input/"
REF=$DIR"Reference/human_g1k_v37.fasta.gz"

function downloadSequence() {
    cd $DIR"Reference/"
    wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
    wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

    gunzip -c $REF >$DIR"Reference/human_g1k_v37.fasta"

    cd $DIR"Input/"
    wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
    wget http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
}

function indexReference() {
    cd $DIR"Reference/"
    $Bwa index $REF
}

function alignReads() {
    $Bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" -t 20 $REF $DIR"Input/ERR194147_1.fastq.gz" $DIR"Input/ERR194147_2.fastq.gz" |
        samtools view -S -b - \
            >$DIR"Input/NA12878.bam"
    $Samtools sort -@ 10 -m 2G $DIR"Input/NA12878.bam" -o $DIR"Input/NA12878_sort.bam"
}

function runFastqc() {
    mkdir -p $DIR"Results/fastqc"
    $Time -v $Fastqc $DIR"Input/ERR194147_1.fastq.gz" $DIR"Input/ERR194147_2.fastq.gz" -o $DIR"Results/fastqc/" -t 6
}

function runBamqc() {
    mkdir -p $DIR"Results/bamqc"
    $Time -v $Bamqc $DIR"Input/NA12878_sort.bam" -o $DIR"Results/bamqc/" -t 6
}

function runPicard() {
    mkdir -p $DIR"Results/picard/"
    $Time -v java -jar $Picard CollectInsertSizeMetrics I=$DIR"Input/NA12878_sort.bam" O=$DIR"Results/picard/out.txt" H=$DIR"Results/picard/out.pdf" VALIDATION_STRINGENCY=SILENT
}

function runHyINDEL() {
    $Time -v $HyINDEL -i $DIR"Input/NA12878_sort.bam" -o $DIR"Output/HyINDEL/" -s 318 -d 78 -l 101 -c 52 -t 22
}

function runLumpy() {
    mkdir -p $DIR"Results/lumpy/"
    $Time -v $Samtools view -@ 4 -b -F 1294 $DIR"Input/NA12878_sort.bam" >$DIR"Input/NA12878_disc.bam"
    $Time -v $Samtools view -@ 4 -h $DIR"Input/NA12878_sort.bam" | extractSplitReads_BwaMem -i stdin | samtools view -Sb - >$DIR"Input/NA12878_split.bam"
    $Time -v lumpyexpress -B $DIR"Input/NA12878_sort.bam" -S $DIR"Input/NA12878_split.bam" -D $DIR"Input/NA12878_disc.bam" -o $DIR"Results/lumpy/NA12878_lumpy.vcf"
}

function runTiddit() {
    mkdir -p $DIR"Results/tiddit/"
    python $Tiddit --sv --bam $DIR"Input/NA12878_sort.bam" -o $DIR"Results/tiddit/NA12878_tiddit"
}

function runSoftSV() {
    mkdir -p $DIR"Results/softsv/"
    $Time -v $SoftSV --noInversions --noDuplications --noTranslocations -s 318 -d 45 -i $DIR"Input/NA12878_sort.bam" -o $DIR"Results/softsv"
}

function runPamir() {
    $Time -v $Pamir --mrsfast-threads 23 --num-worker 23 --donot-remove-contaminant \
        -p $DIR"Results/pamir_aln" -r $DIR"Reference/human_g1k_v37.fasta" \
        --files alignment=$DIR"Input/NA12878_readnamesorted.bam"
}

function runPopins() {
    mkdir -p $DIR"Output/popins" && cd $DIR"Output/popins"
    ln -s $DIR"Reference/human_g1k_v37.fasta" genome.fa
    ln -s $DIR"Reference/human_g1k_v37.fasta.fai" genome.fa.fai
    $Popins assemble --sample sample1 $DIR"Input/NA12878_sort.bam" -t 22 -m 2G
    $Popins merge
    $Popins contigmap sample1 -t 22 -m 2G
    $Popins place-refalign
    $Popins place-splitalign sample1
    $Popins place-finish
    $Popins genotype sample1
}

# downloadSequence
# indexReference
# alignReads

# runFastqc
# runBamqc

# runPicard

# runLumpy
# runTiddit
# runPamir
# runSoftSV
runHyINDEL
# runPopins
