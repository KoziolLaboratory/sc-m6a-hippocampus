#!/bin/bash
#SBATCH -J multiqc
#SBATCH -o multiqc_MeRIP_%A.out
#SBATCH -e multiqc_MeRIP_%A.err
#SBATCH -p q_cn
#SBATCH -n 1

module load multiqc

# Parameters of the current job
pwd; hostname; date

# Input arguments
datadir=$1

multiqc $datadir

