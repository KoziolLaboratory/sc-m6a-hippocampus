#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=align_%A.out
#SBATCH --error=align_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=160G

# Description:  Here I am aligning bulk dart-seq mice hippocampus reads to mm39 reference genome using bwa-mem

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
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf mm39.chromFa.tar.gz
cat chr* > mm39.fa
rm chr*
rm mm39.chromFa.tar.gz
cp /home/Magda_lab/maitenat/DATA/indeces/mm* ${SCRATCH_DIRECTORY}

# Align the data
bwa mem mm39.fa ${outdir}/${samplename}_joint.fq.gz | samtools sort -O BAM -o ${outdir}/${samplename}_joint.bam -

# After the job is done we copy our output back to our directory
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
rm ${outdir}/${samplename}_joint.fq.gz
