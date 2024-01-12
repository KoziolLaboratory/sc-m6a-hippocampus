#!/bin/bash
#SBATCH -J cutadapt
#SBATCH -o cutadapt_%A.out
#SBATCH -e cutadapt_%A.err
#SBATCH -p q_cn
#SBATCH -n 1

module load cutadapt

# Parameters of the current job
pwd; hostname; date

# Adapter sequences
adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Input fq files
input_r1=$1
input_r2=$2

output_r1=${input_r1/.fq.gz/_clean.fq.gz}
output_r2=${input_r2/.fq.gz/_clean.fq.gz}

# Run cutadapt to remove adapters and write to temporary files
cutadapt -a $adapter1 -A $adapter2 -m 15 -j 3 -p temp_r2 -o temp_r1 $input_r1 $input_r2

# Run cutadapt again to trim 10 bases from the start of each read
cutadapt -u 10 -U 10 -o $output_r1 -p $output_r2 temp_r1 temp_r2

# Remove the temporary files
rm temp_r1 temp_r2
