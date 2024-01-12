#!/bin/bash
#SBATCH -J peak_intersection
#SBATCH -o peak_intersection_%A.out
#SBATCH -p q_cn
#SBATCH -n 1

module load bedtools/2.29.2

# Parameters of the current job
pwd; hostname; date

# (Narrow) peak files
peak1=$1
peak2=$2
peak3=$3
outdir=`dirname $peak1`

# Individual intersection
bedtools intersect -a $peak1 -b $peak2 -f 0.5 -F 0.5 -e > rep1_rep2_intersect.bed
bedtools intersect -a $peak1 -b $peak3 -f 0.5 -F 0.5 -e > rep1_rep3_intersect.bed
bedtools intersect -a $peak2 -b $peak3 -f 0.5 -F 0.5 -e > rep2_rep3_intersect.bed

# Intersection of all three
bedtools intersect -a $peak1 -b $peak2 -f 0.5 -F 0.5 -e | bedtools intersect -a stdin -b $peak3 -f 0.5 -F 0.5 -e > reps_intersect.bed

# Copy results to output directory
cp rep1_rep2_intersect.bed $outdir
cp rep1_rep3_intersect.bed $outdir
cp rep2_rep3_intersect.bed $outdir
cp reps_intersect.bed $outdir