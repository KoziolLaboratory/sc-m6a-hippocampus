#!/bin/bash

# In this script we are annotating the genes to the MeRIP sites

module load bedtools/2.29.2

# Peak files of each replicate
merip_rep1=$1
merip_rep2=$2
merip_rep3=$3
# Replicate intersection bed file
yth_merip=$4
# GTF file
gtf=/GPFS/Magda_lab_permanent/maitenat/mm39.ncbiRefSeq.gtf

# Annotate the gene symbol to the MeRIP bed file using gtf file
# First, we need to convert the gtf file to bed file
# Then, we need to intersect the MeRIP bed file with the gtf bed file
# Then, we need to remove the duplicates

# Convert gtf to bed file
awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$7"\t"$9}' $gtf > mm39.ncbiRefSeq.bed

# Intersect MeRIP bed file with gtf bed file
bedtools intersect -a $yth_merip -b mm39.ncbiRefSeq.bed -wa -wb > intersect_merip_annot.bed

# Remove duplicates based on the 3 first columns (pos) and gene (14th column)
awk '!seen[$1,$2,$3,$14]++' intersect_merip_annot.bed > intersect_merip_annot_nodup.bed

# Not all regions in the MeRIP have been annotated!
# Get difference between non-annotated and annotated files
# This will give us the non-annotated genes
awk 'NR==FNR{a[$1,$2,$3];next} !($1,$2,$3) in a' intersect_merip_annot_nodup.bed $yth_merip > intersect_merip_nonannot.bed

# We will also annotate each of the replicates to make venn diagrams with the genes afterwards
# Intersect MeRIP rep1 file with gtf bed file
bedtools intersect -a $merip_rep1 -b mm39.ncbiRefSeq.bed -wa -wb > merip_rep1_annot.bed
# Remove duplicates based on the 3 first columns (pos) and gene (14th column)
awk '!seen[$1,$2,$3,$14]++' merip_rep1_annot.bed > merip_rep1_annot_nodup.bed
# Intersect MeRIP rep2 file with gtf bed file
bedtools intersect -a $merip_rep2 -b mm39.ncbiRefSeq.bed -wa -wb > merip_rep2_annot.bed
# Remove duplicates based on the 3 first columns (pos) and gene (14th column)
awk '!seen[$1,$2,$3,$14]++' merip_rep2_annot.bed > merip_rep2_annot_nodup.bed
# Intersect MeRIP rep3 file with gtf bed file
bedtools intersect -a $merip_rep3 -b mm39.ncbiRefSeq.bed -wa -wb > merip_rep3_annot.bed
# Remove duplicates based on the 3 first columns (pos) and gene (14th column)
awk '!seen[$1,$2,$3,$14]++' merip_rep3_annot.bed > merip_rep3_annot_nodup.bed
