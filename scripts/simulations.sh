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

THREADS=22

Time="/usr/bin/time"
SVsim="/home/iiit/alok.t/Tools/SVsim/SVsim"
Art="/home/iiit/alok.t/Tools/art_src_MountRainier_Linux/art_illumina"
Bwa="/home/iiit/alok.t/miniconda2/bin/bwa"
Samtools="/home/iiit/alok.t/miniconda2/bin/samtools"
Fastqc="/home/iiit/alok.t/Tools/FastQC/fastqc"
Bamqc="/home/iiit/alok.t/Tools/BamQC/bin/bamqc"
Picard="/home/iiit/alok.t/Tools/picard.jar"
Lumpy="/home/iiit/alok.t/miniconda2/bin/lumpyexpress"
Lumpy_extract="/home/iiit/alok.t/miniconda2/bin/extractSplitReads_BwaMem"
Tiddit="/home/iiit/alok.t/Tools/TIDDIT-TIDDIT-2.6.0/TIDDIT.py"
SoftSV="/home/iiit/alok.t/Tools/SoftSV_1.4.2/SoftSV"
Popins="/home/iiit/alok.t/Tools/popins/popins"
Pamir="/home/iiit/alok.t/Tools/pamir-master/pamir.py"
HyINDEL_path="/home/iiit/alok.t/Tools/HyINDEL/scripts/"
HyINDEL="/home/iiit/alok.t/install/HyINDEL/bin/HyINDEL"

DIR="/scratch/alok/Simulations/"
mkdir -p $DIR"Input/" $DIR"Reference/" $DIR"Output/"
REF=$DIR"Reference/human_g1k_v37.fasta"

function downloadReference() {
    wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz \
        -O $DIR"Reference/human_g1k_v37.fasta.gz"
    gunzip -c $DIR"Reference/human_g1k_v37.fasta.gz"
}

function indexReference() {
    $Samtools faidx $REF
    $Bwa index $REF
}

function generateSample() {
    # Generate variants for tmp sample
    python $HyINDEL_path/helperSimulations.py genPositions >$DIR"Input/tmp.txt"
    # Insert variants in reference to create a tmp sample
    python $SVsim -i $DIR"Input/tmp.txt" -r $REF -o $DIR"Input/tmp" -W -d -s 572239

    # Read positions and generate variants for final sample
    python $HyINDEL_path/helperSimulations.py genVariants $DIR"Input/"

    # Insert variants into sample1
    python $SVsim -i $DIR"Input/1.txt" -r $REF -o $DIR"Input/1" -W -d -s 918487
    # Insert variants into sample2
    python $SVsim -i $DIR"Input/2.txt" -r $REF -o $DIR"Input/2" -W -d -s 749933
    # Merge first, second sample to create diploid sample
    cat $DIR"Input/1.fasta" $DIR"Input/2.fasta" >$DIR"Input/sim.fasta"
    $Samtools faidx $DIR"Input/sim.fasta"
}

# Args: Cov, Cov/2
function generateAlignment() {
    ALN_FOLDER=$DIR"Input/"$1"x/"
    mkdir -p $ALN_FOLDER

    # Generate reads
    $Time -v $Art -ss HS25 -na -q -i $DIR"Input/sim.fasta" -d $1"x_sim" -p -l 100 -f $2 -m 350 -s 20 -o $ALN_FOLDER"sim"

    # Align reads
    $Time -v $Bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" -t $THREADS $REF $ALN_FOLDER"sim1.fq" $ALN_FOLDER"sim2.fq" >$ALN_FOLDER"sim.sam"

    # Sort and convert to bam
    $Time -v $Samtools sort -@ $THREADS -m 2G $ALN_FOLDER"sim.sam" -o $ALN_FOLDER"sim.bam"
    # Index bam
    $Time -v $Samtools index $ALN_FOLDER"sim.bam"
}

# Args: Cov
function runPicard() {
    mkdir -p $DIR"Output/"$1"x/picard"
    ALN_FOLDER=$DIR"Input/"$1"x/"
    $Time -v java -jar $Picard CollectInsertSizeMetrics I=$ALN_FOLDER"sim.bam" O=$DIR"Output/"$1"x/picard/out.txt" H=$DIR"Output/"$1"x/picard/out.pdf" VALIDATION_STRINGENCY=SILENT
}

# Args: Cov
function runLumpy() {
    mkdir -p $DIR"Output/"$1"x/lumpy"
    ALN_FOLDER=$DIR"Input/"$1"x/"
    $Samtools view -b -F 1294 -@ 4 $ALN_FOLDER"sim.bam" >$ALN_FOLDER"sim_disc.bam"
    $Samtools view -h -@ 4 $ALN_FOLDER"sim.bam" | $Lumpy_extract -i stdin | samtools view -Sb - >$ALN_FOLDER"sim_split.bam"
    $Lumpy -B $ALN_FOLDER"sim.bam" -S $ALN_FOLDER"sim_split.bam" -D $ALN_FOLDER"sim_disc.bam" -o $DIR"Output/"$1"x/lumpy/lumpy_sim_"$1"x.vcf"
}

# Args: Cov
function runTiddit() {
    mkdir -p $DIR"Output/"$1"x/tiddit"
    ALN_FOLDER=$DIR"Input/"$1"x/"
    python $Tiddit --sv --bam $ALN_FOLDER"sim.bam" -o $DIR"Output/"$1"x/tiddit/tiddit_sim_"$1"x"
}

# Args: Cov
function runSoftSV() {
    mkdir -p $DIR"Output/"$1"x/softsv"
    ALN_FOLDER=$DIR"Input/"$1"x/"
    $SoftSV -i $ALN_FOLDER"sim.bam" -o $DIR"Output/"$1"x/softsv"
}

# Args: Cov
function runPopins() {
    mkdir -p $DIR"Output/"$1"x/popins" && cd $DIR"Output/"$1"x/popins"
    ALN_FOLDER=$DIR"Input/"$1"x/"
    ln -s $DIR"Reference/human_g1k_v37.fasta" genome.fa
    ln -s $DIR"Reference/human_g1k_v37.fasta.fai" genome.fa.fai
    $Popins assemble --sample sample1 $ALN_FOLDER"sim.bam" -t $THREADS -m 2G
    $Popins merge
    $Popins contigmap sample1 -t $THREADS -m 2G
    $Popins place-refalign
    $Popins place-splitalign sample1
    $Popins place-finish
    $Popins genotype sample1
}

# Args: Cov
function runPamir() {
    ALN_FOLDER=$DIR"Input/"$1"x/"
    $Pamir --mrsfast-threads $THREADS --num-worker $THREADS --donot-remove-contaminant \
        -p $DIR"Output/"$1"x/pamir_mrsfast_best" -r $REF \
        --files mrsfast-best-search=$ALN_FOLDER"sim.bam"
}

# Args: Cov
function runHyINDEL() {
    ALN_FOLDER=$DIR"Input/"$1"x/"
    $Time -v $HyINDEL -i $ALN_FOLDER"sim.bam" -o $DIR"Output/"$1"x/hyindel/" -s 350 -d 20 -l 100 -c $1 -t $THREADS
}

# downloadReference
# indexReference
# generateSample
# generateAlignment 10 5
# generateAlignment 20 10
# generateAlignment 30 15

# runPicard 30
# runPicard 20
# runPicard 10

# runHyINDEL 30
# runHyINDEL 20
# runHyINDEL 10

# runLumpy 10
# runLumpy 20
# runLumpy 30

# runTiddit 10
# runTiddit 20
# runTiddit 30

# runSoftSV 10
# runSoftSV 20
# runSoftSV 30

# runPopins 10
# runPopins 20
# runPopins 30

# runPamir 30
# runPamir 20
# runPamir 10
