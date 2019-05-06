#!/bin/bash
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p long
#SBATCH -w node56
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-type=END

# Add / at end for directory path
SIM_DIR="/scratch/alok/Simulations/"
REF=$SIM_DIR"Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa"

# Tools
Bwa="/home/iiit/alok.t/miniconda2/bin/bwa"
Samtools="/home/iiit/alok.t/miniconda2/bin/samtools"
# https://github.com/GregoryFaust/SVsim
SVsim="/home/iiit/alok.t/Tools/SVsim/SVsim"
# https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
Art="/home/iiit/alok.t/Tools/art_src_MountRainier_Linux/art_illumina"

Lumpy="/home/iiit/alok.t/miniconda2/bin/lumpyexpress"
Lumpy_extract="/home/iiit/alok.t/miniconda2/bin/extractSplitReads_BwaMem"
Tiddit="/home/iiit/alok.t/Tools/TIDDIT-TIDDIT-2.6.0/TIDDIT.py"
SoftSV="/home/iiit/alok.t/Tools/SoftSV_1.4.2/SoftSV"
PlusMinus="/home/iiit/alok.t/install/indel-detect/bin/plusminus"

mkdir -p $SIM_DIR"Reference"
mkdir -p $SIM_DIR"Results"

function downloadSequence() {
    wget http://ftp.ensembl.org/pub/grch37/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz -O $SIM_DIR"Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz"
    gunzip -c $SIM_DIR"Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz" >$REF
    $Samtools faidx $REF
    $Bwa index $REF

    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed -O $SIM_DIR"Reference/Personalis_1000_Genomes_deduplicated_deletions.bed"
}

function generateSample() {
    python $(dirname $0)/range.py $SIM_DIR"Reference/Personalis_1000_Genomes_deduplicated_deletions.bed" >$SIM_DIR"mod.txt"

    python $SVsim -i $SIM_DIR"mod.txt" -r $REF -o $SIM_DIR"mod" -W -d
    rm $SIM_DIR"mod.fasta"

    python $(dirname $0)/half.py $SIM_DIR

    # Insert variants into reference
    python $SVsim -i $SIM_DIR"sim1.txt" -r $REF -o $SIM_DIR"sim1" -W -d
    python $SVsim -i $SIM_DIR"sim2.txt" -r $REF -o $SIM_DIR"sim2" -W -d
}

function generateAlignment() {
    mkdir -p $SIM_DIR$1"x/"

    # Generate reads
    $Art -ss HS20 -na -q -i $SIM_DIR"sim1.fasta" -d sim1 -p -l 100 -f $2 -m 350 -s 20 -o $SIM_DIR$1"x/sim1_"
    $Art -ss HS20 -na -q -i $SIM_DIR"sim2.fasta" -d sim2 -p -l 100 -f $2 -m 350 -s 20 -o $SIM_DIR$1"x/sim2_"

    # Merge reads
    cat $SIM_DIR$1"x/sim2_1.fq" >>$SIM_DIR$1"x/sim1_1.fq"
    cat $SIM_DIR$1"x/sim2_2.fq" >>$SIM_DIR$1"x/sim1_2.fq"

    # Align reads
    $Bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" -t 22 $REF $SIM_DIR$1"x/sim1_1.fq" $SIM_DIR$1"x/sim1_2.fq" >$SIM_DIR$1"x/sim.sam"

    # Sort and convert to bam
    $Samtools sort -@ 10 -M 2G $SIM_DIR$1"x/sim.sam" -o $SIM_DIR$1"x/sim.bam"
    # Index bam
    $Samtools index $SIM_DIR$1"x/sim.bam"

    # Remove files
    rm $SIM_DIR$1"x/sim2_1.fq" $SIM_DIR$1"x/sim2_2.fq" $SIM_DIR$1"x/sim.sam"
}

function runLumpy() {
    mkdir -p $SIM_DIR"Results/"$1"x/lumpy"
    $Samtools view -b -F 1294 -@ 4 $SIM_DIR$1"x/sim.bam" >$SIM_DIR$1"x/sim_disc.bam"
    $Samtools view -h -@ 4 $SIM_DIR$1"x/sim.bam" | $Lumpy_extract -i stdin | samtools view -Sb - >/$SIM_DIR$1"x/sim_split.bam"
    $Lumpy -B $SIM_DIR$1"x/sim.bam" -S $SIM_DIR$1"x/sim_split.bam" -D $SIM_DIR$1"x/sim_disc.bam" -o $SIM_DIR"Results/"$1"x/lumpy/lumpy_sim_"$1"x.vcf"
}

function runTiddit() {
    mkdir -p $SIM_DIR"Results/"$1"x/tiddit"
    python $Tiddit --sv --bam $SIM_DIR$1"x/sim.bam" -o $SIM_DIR"Results/"$1"x/tiddit/tiddit_sim_"$1"x" -i 410
}

function runSoftSV() {
    mkdir -p $SIM_DIR"Results/"$1"x/softsv"
    $SoftSV -i $SIM_DIR$1"x/sim.bam" -o $SIM_DIR"Results/"$1"x/softsv"
}

function runPlusMinus() {
    mkdir -p $SIM_DIR"Results/"$1"x/plusMinus"
    $PlusMinus -i $SIM_DIR$1"x/sim.bam" -o $SIM_DIR"Results/"$1"x/plusMinus" -s 350 -d 20 -l 100 -c $1 -t 8
}

# Uncomment lines below to run

# downloadSequence

# generateSample

# Arguments: Cov, Cov/2
# generateAlignment 5 2.5
# generateAlignment 10 5
# generateAlignment 20 10
# generateAlignment 30 15

# runLumpy 5
# runLumpy 10
# runLumpy 20
# runLumpy 30

runTiddit 5
runTiddit 10
runTiddit 20
runTiddit 30

# runSoftSV 5
# runSoftSV 10
# runSoftSV 20
# runSoftSV 30

# runPlusMinus 5
# runPlusMinus 10
# runPlusMinus 20
# runPlusMinus 30
