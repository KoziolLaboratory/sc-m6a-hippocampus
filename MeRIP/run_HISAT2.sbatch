#!/bin/bash
#SBATCH -J alignment_byHISAT2
#SBATCH -o alighment_byHISAT2.out
#SBATCH -p q_cn
#SBATCH -n 1

module load samtools
module load hisat/2.1.0

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
fq1=$1
fq2=$2
outdir=/home/Magda_lab/maitenat/scratch60/MeRIP/aligned

# Copy necessary files to scratch directory
cp $fq1 ${SCRATCH_DIRECTORY}
cp $fq2 ${SCRATCH_DIRECTORY}
fq1_basename=`basename $fq1`
fq2_basename=`basename $fq2`
samplename="${fq1_basename%%.*}"

# Build index
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf mm39.chromFa.tar.gz
cat chr* > mm39.fa
rm chr*
rm mm39.chromFa.tar.gz
hisat2-build mm39.fa mm39

hisat2 -p 10 -x mm39 -1 ${fq1_basename} -2 ${fq2_basename} 2>${samplename}.log | samtools view -@ 16 -bS -f 2 -F 256 -q 30 | samtools sort -@ 10 -o ${samplename}.Hisat_aln.sorted.bam
samtools index ${samplename}.Hisat_aln.sorted.bam

# After the job is done we copy our output back to our directory
cp ${samplename}.Hisat_aln.sorted.bam ${outdir}
cp ${samplename}.Hisat_aln.sorted.bam.bai ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}


