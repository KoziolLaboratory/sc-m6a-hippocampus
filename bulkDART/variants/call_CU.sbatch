#!/bin/bash
#SBATCH --job-name=call_CU
#SBATCH --output=call_CU_%A.out
#SBATCH --error=call_CU_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_fat_l
#SBATCH --mem-per-cpu=140G

# Description:  Here I am using CTK to call C>U deaminations. Before that, I remove secondary alignments and mark duplicates. 
# The final C>U variant set is called after merging 3 replicate sets. 

module load ctk/1.1.4
module load picard/2.20.4
module load samtools/1.12
module load python/3.7.6

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
bamfile=$1

# Copy necessary files to scratch directory
cp ${bamfile} ${SCRATCH_DIRECTORY}
bamfile_basename=`basename ${bamfile}`
samplename=`echo "${bamfile_basename//_1_joint.bam/}"`

# Remove secondary alignments
samtools view -b -F 2048,256 ${bamfile_basename} > ${samplename}_nosec_markdup.bam

# Remove duplicates
picard MarkDuplicates \
REMOVE_DUPLICATES=true \
I=${samplename}_nosec_markdup.bam \
O=${samplename}_markdup.bam \
M=marked_dup_metrics.txt

rm ${samplename}_nosec_markdup.bam

# Identify C-U mutations using CIMS
# First need to convert bam to sam
samtools view -h -o ${samplename}_markdup.sam ${samplename}_markdup.bam
# Convert bam to bed
parseAlignment.pl -v --map-qual 1 --min-len 18 --mutation-file ${samplename}.mutation.txt ${samplename}_markdup.sam ${samplename}_markdup.tag.bed

# Collapse PCR duplicates with CTK 
tag2collapse.pl -v -big -weight --weight-in-name --keep-max-score --keep-tag-name -c ${SCRATCH_DIRECTORY}/cache ${samplename}_markdup.tag.bed ${samplename}_markdup.tag.uniq.bed

# Get mutations in unique tags
python /usr/nzx-cluster/apps/ctk/ctk/joinWrapper.py ${samplename}.mutation.txt ${samplename}_markdup.tag.uniq.bed 4 4 N ${samplename}.tag.uniq.mutation.txt

# Run CIMS
# Filter for C>T (G>A) substitutions
awk '{if($6=="+" && $8=="C" && $9==">" && $10=="T" || $6=="-" && $8=="G" && $9==">" && $10=="A") {print $0}}' ${samplename}.tag.uniq.mutation.txt | cut -f 1-6 > ${samplename}.tag.uniq.subs.bed

# Run CIMS
CACHEHOME=${SCRATCH_DIRECTORY}
export CACHEHOME
CIMS.pl -big -v -n 5 -p -c ${SCRATCH_DIRECTORY}/cache ${samplename}_markdup.tag.uniq.bed ${samplename}.tag.uniq.subs.bed ${samplename}.cims

# After the job is done we copy our output back to our directory
outdir=/home/Magda_lab/maitenat/scratch60/dartseq_replicates/aligned
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}.cims ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
