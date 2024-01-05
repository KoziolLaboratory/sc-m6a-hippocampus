#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --output=rename_%A.out
#SBATCH --error=rename_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=20G

# Description:  Here I am renaming the fastq files to match the format of the bam files

# Parameters of the current job
pwd; hostname; date

outdir=/home/Magda_lab/maitenat/scratch60/dartseq_replicates/
fq_file=$1
fq_basename=`basename ${fq_file}`
samplename="${fq_basename%%.*}"

# for read 2, change :1 to :2
zcat ${fq_file} | awk '{if (NR % 4 == 1) {print $1":1 " $2} else {print $0} }' > ${outdir}/${samplename}_rename.fastq
gzip ${outdir}/${samplename}_rename.fastq
rm ${outdir}/${samplename}_rename.fastq