#!/bin/bash
#SBATCH --job-name=filtbull
#SBATCH --output=filter_bullseye_%A.out
#SBATCH --error=filter_bullseye_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=80G

module load bedtools/2.29.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
sites_bed=$1
barcode_list=$2
snp_bedfile=mgp_REL2021_snps.bed
eva_snp_bedfile=GCA_000001635.9_current_ids.bed
apobec_bulk_23=B26-B27-B28_2out3.bed
wt_bulk_23=B17-B18-B19_2out3.bed

# Copy necessary files to scratch directory
cp $sites_bed ${SCRATCH_DIRECTORY}
bed_basename=`basename ${sites_bed}`
cp $barcode_list ${SCRATCH_DIRECTORY}
barcode_basename=`basename ${barcode_list}`
cp $snp_bedfile ${SCRATCH_DIRECTORY}
snp_basename=`basename ${snp_bedfile}`
cp $eva_snp_bedfile ${SCRATCH_DIRECTORY}
eva_snp_basename=`basename ${eva_snp_bedfile}`
outdir=`dirname $sites_bed`
samplename="${bed_basename%%.*}"

# Discard those variants with 0 barcode
awk '$12 != 0 {print $0}' $bed_basename > filt1.bed

# Filter snps
bedtools subtract -A -a filt1.bed -b ${snp_basename} > filt2.bed
bedtools subtract -A -a filt2.bed -b ${eva_snp_basename} > filt3.bed

# Remove from sc sites that are in (bulk) apobec and wt
bedtools intersect -a filt3.bed -b $apobec_bulk_23 $wt_bulk_23 -v > filt4.bed

n_filt1=`wc -l filt4.bed | awk '{print $1}'`
echo 'There are' $((n_filt1)) 'total variants discarding SNPs and variants in apobec and wt.'

# Keep only those variants whose barcode is in a tracked cell
awk 'NR==FNR{A[$1];next}($12 in A)' ${barcode_basename} filt4.bed > ${samplename}_filtered_all.bed
n_filt2=`wc -l ${samplename}_filtered_all.bed | awk '{print $1}'`
echo 'Out of the previous, there are' $((n_filt2)) 'total variants assigned to a cell.'

# After the job is done we copy our output back to our directory
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}_filtered_all.bed ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
