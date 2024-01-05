#!/bin/bash
#SBATCH --job-name=filtCIMS
#SBATCH --output=filter_CIMS_%A.out
#SBATCH --error=filter_CIMS_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=40G

# Description:  Here we are filtering the results of running CIMS on bulk dart-seq mice hippocampus reads
# #

module load bedtools/2.29.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
cims_bedfile=$1

# Copy necessary files to scratch directory
cp $cims_bedfile ${SCRATCH_DIRECTORY}
cims_bed_basename=`basename ${cims_bedfile}`

snp_bedfile=/GPFS/Magda_lab_permanent/maitenat/mgp_REL2021_snps.bed
eva_snp_bedfile=/GPFS/Magda_lab_permanent/maitenat/GCA_000001635.9_current_ids.bed
cp $snp_bedfile ${SCRATCH_DIRECTORY}
snp_basename=`basename ${snp_bedfile}`
cp $eva_snp_bedfile ${SCRATCH_DIRECTORY}
eva_snp_basename=`basename ${eva_snp_bedfile}`

samplename="${cims_bed_basename%%.*}"

# Filter by m >= 2, fdr > 1, k>= 10 
awk "{if(\$9<1 && \$8>=2 && \$7>=10) {print \$0}}" ${cims_bed_basename} > filt0.cims
# Add m/k
awk -v OFS='\t' '{a=$8/$7;print $0,a;}' filt0.cims > filt1.cims
# Filter snps
bedtools subtract -A -a filt1.cims -b ${snp_basename} > filt2.cims
bedtools subtract -A -a filt2.cims -b ${eva_snp_basename} > filt3.cims

# Filter by m/k
awk "{if(\$11>=0.05) {print \$0}}" filt3.cims > "${samplename}_filtered.cims"

# After the job is done we copy our output back to our directory
outdir=/GPFS/Magda_lab_temp/maitenat/dartseq_replicates/aligned
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/*filtered.cims ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
