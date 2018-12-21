for c in `seq 1 22` X Y
do

# large
sort -k1,1 -k2,2n large/${c}_filter_dels_large.txt > ${c}_l0.bed
awk -f filter_support_large.awk ${c}_l0.bed > ${c}_l1.bed
samtools bedcov ${c}_l1.bed /Volumes/HDD/Data/30x/chr/${c}.bam > ${c}_l2.bed
awk -f filter_coverage.awk ${c}_l2.bed > ${c}_l3.bed
bedtools merge -i ${c}_l3.bed > ${c}_large.bed

# small
sort -k1,1 -k2,2n small/${c}_filter_dels_small.txt > ${c}_s0.bed
awk -f filter_support_small.awk ${c}_s0.bed > ${c}_s1.bed
samtools bedcov ${c}_s1.bed /Volumes/HDD/Data/30x/chr/${c}.bam > ${c}_s2.bed
awk -f filter_coverage.awk ${c}_s2.bed > ${c}_s3.bed
bedtools merge -i ${c}_s3.bed > ${c}_small.bed

# split
sort -k1,1 -k2,2n split/${c}_filter_deletions.txt > ${c}_sr0.bed
awk -f filter_support_split.awk ${c}_sr0.bed > ${c}_sr1.bed
bedtools merge -i ${c}_sr1.bed > ${c}_sr2.bed
samtools bedcov ${c}_sr2.bed /Volumes/HDD/Data/30x/chr/${c}.bam > ${c}_sr3.bed
awk -f filter_coverage.awk ${c}_sr3.bed > ${c}_split.bed


cat ${c}_large.bed ${c}_small.bed > ${c}.bed
rm -f ${c}_l0.bed ${c}_l1.bed ${c}_l2.bed ${c}_l3.bed ${c}_large.bed ${c}_s0.bed ${c}_s1.bed ${c}_s2.bed ${c}_s3.bed ${c}_small.bed ${c}_sr0.bed ${c}_sr1.bed ${c}_sr2.bed ${c}_sr3.bed

echo ${c} "done"
done