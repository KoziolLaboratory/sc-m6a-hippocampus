#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=align_%A.out
#SBATCH --error=align_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=160G

# Description:  Here I am aligning bulk dart-seq HEK cell-line reads to hg19 reference genome using bwa-mem

module load bwa/0.7.17
module load samtools/1.12

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
fq1=$1
fq2=$2
outdir=/home/Magda_lab/maitenat/scratch60/dartseq_replicates/aligned

# Copy necessary files to scratch directory
fq1_basename=`basename $fq1`
samplename="${fq1_basename%%.*}"
cat $fq1 $fq2 > ${outdir}/${samplename}_joint.fq.gz
cp /home/Magda_lab/maitenat/DATA/chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf chromFa.tar.gz
cat chr* > hg19.fa
rm chr*
rm chromFa.tar.gz
cp /home/Magda_lab/maitenat/DATA/indeces/hg* ${SCRATCH_DIRECTORY}

# Align the data
bwa mem hg19.fa ${outdir}/${samplename}_joint.fq.gz | samtools sort -O BAM -o ${outdir}/${samplename}_joint.bam -

# After the job is done we copy our output back to our directory
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
rm ${outdir}/${samplename}_joint.fq.gz
