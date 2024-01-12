#!/bin/bash
#SBATCH -J fastqc
#SBATCH -o fastqc_MeRIP_%A.out
#SBATCH -e fastqc_MeRIP_%A.err
#SBATCH -p q_cn
#SBATCH -n 1

module load fastqc

# Parameters of the current job
pwd; hostname; date

# Input arguments
R1=$1
R2=$2
outdir=$3

fastqc -t 12 -o $outdir $R1 $R2
