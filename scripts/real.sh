#!/bin/bash
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p long
#SBATCH -w node35
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --mail-type=END

function runLumpy() {
    mkdir -p /scratch/alok/Real/Results/lumpy
    samtools view -b -F 1294 -@ 4 /scratch/alok/Real/RMNISTHS_30xdownsample.bam >/scratch/alok/Real/30x_disc.bam
    samtools view -h -@ 4 /scratch/alok/Real/RMNISTHS_30xdownsample.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb - >/scratch/alok/Real/30x_split.bam
    lumpyexpress -B /scratch/alok/Real/RMNISTHS_30xdownsample.bam -S /scratch/alok/Real/30x_split.bam -D /scratch/alok/Real/30x_disc.bam -o /scratch/alok/Real/Results/lumpy/lumpy_real.vcf
}

function runTiddit() {
    mkdir -p /scratch/alok/Real/Results/tiddit
    python /home/iiit/alok.t/Tools/TIDDIT-TIDDIT-2.6.0/TIDDIT.py --sv --bam /scratch/alok/Real/RMNISTHS_30xdownsample.bam -o /scratch/alok/Real/Results/tiddit/tiddit_real
}

function runSoftSV() {
    mkdir -p /scratch/alok/Real/Results/softsv
    /home/iiit/alok.t/Tools/SoftSV_1.4.2/SoftSV --noInversions --noDuplications --noTranslocations -s 550 -d 100 -i /scratch/alok/Real/RMNISTHS_30xdownsample.bam -o /scratch/alok/Real/Results/softsv
}

function runPlusMinus() {
    mkdir -p /scratch/alok/Real/Results/plusminus
    /home/iiit/alok.t/install/indel-detect/bin/plusminus -i /scratch/alok/Real/RMNISTHS_30xdownsample.bam -o /scratch/alok/Real/Results/plusminus -s 550 -d 100 -l 148 -c 33 -t 22
}

# runLumpy
runTiddit
# runSoftSV
# runPlusMinus
