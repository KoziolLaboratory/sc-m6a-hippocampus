#!/bin/bash
#SBATCH -J rRNA_remove
#SBATCH -o rRNA_remove_%A.out
#SBATCH -e rRNA_remove_%A.err
#SBATCH -p q_cn
#SBATCH -n 1

module load bowtie/2.3.5

fq1=$1
fq2=$2
fq1_basename=`basename $fq1`
samplename="${fq1_basename%%.*}"
outname=${samplename}_rRNA_removed

rRNA_seq=/home/Magda_lab/maitenat/DATA/GCF_000001635.27.p13_rna_from_genomic.fna 

# Obtain the id of rRNA sequence
cat GCF_000001405.25_GRCh37.p13_rna_from_genomic.fna | grep "^>" | grep "gbkey=rRNA" | awk '{print $1}' | sed 's/>//g' > id.list

# Obtain the sequence of rRNA sequence
seqkit grep -f id.list GCF_000001405.25_GRCh37.p13_rna_from_genomic.fna > rRNA.fa

# Build index of rRNA.fa
bowtie2-build rRNA.fa rRNA

# Remove rRNA using bowtie2
bowtie2 --very-sensitive-local --no-unal -I 1 -X 100 -p 6 -x rNA -1 $fq1 -2 $fq2 --un-conc-gz $outname
