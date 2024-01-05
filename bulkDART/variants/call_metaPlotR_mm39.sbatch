#!/bin/bash
#SBATCH --job-name=metaplot
#SBATCH --output=metaplotr_%A.out
#SBATCH --error=metaplotr_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=80G

# Description:  Here we are calling metaPlotR on mice hippocampus CIMS data
# #

module load perl/5.34.0
module load bedtools/2.29.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
mkdir -p genome
cims_bedfile=$1

# Copy necessary files to scratch directory
cp $cims_bedfile ${SCRATCH_DIRECTORY}
cims_bed_basename=`basename $cims_bedfile`
samplename="${cims_bed_basename%%.*}"
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}/genome
cd genome
tar -xzvf mm39.chromFa.tar.gz
cd ..
cp /home/Magda_lab/maitenat/code/metaPlotR-master.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf metaPlotR-master.tar.gz
cp /GPFS/Magda_lab_permanent/maitenat/GRCm39_gencode.genePred ${SCRATCH_DIRECTORY}

# Sort bedfile
# Select just first 6 columns of bedfile
cut -f 1-6 ${cims_bed_basename} > file.bed 
sort -k1,1 -k2,2n file.bed > ${samplename}.sorted.bed

# Annotate reference genome
perl metaPlotR-master/make_annot_bed.pl --genomeDir genome/ --genePred GRCm39_gencode.genePred > GRCm39_annot.bed
sort -k1,1 -k2,2n GRCm39_annot.bed > GRCm39_sorted.annot.bed

# Catalogue transcriptome
perl metaPlotR-master/size_of_cds_utrs.pl --annot GRCm39_sorted.annot.bed > region_sizes.txt

# Annotate CIMS bedfile
perl metaPlotR-master/annotate_bed_file.pl --bed ${samplename}.sorted.bed --bed2 GRCm39_sorted.annot.bed > annot_${samplename}.bed

# Metagene coordinates
perl metaPlotR-master/rel_and_abs_dist_calc.pl --bed annot_${samplename}.bed --regions region_sizes.txt > ${samplename}.dist.measures.txt

# After the job is done we copy our output back to our directory
outdir=/GPFS/Magda_lab_temp/maitenat/dartseq_replicates/aligned/metaplotR
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}.dist.measures.txt ${outdir}
cp -r ${SCRATCH_DIRECTORY}/annot_${samplename}.bed ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
