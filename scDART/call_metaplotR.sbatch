#!/bin/bash
#SBATCH --job-name=metaplotr
#SBATCH --output=metaplotr_%A.out
#SBATCH --error=metaplotr_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=80G

# Description:  Here we are calling metaPlotR on all sites

module load perl/5.34.0
module load bedtools/2.29.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
mkdir -p genome
bedfile=$1
outdir=$2

# Copy necessary files to scratch directory
cp $bedfile ${SCRATCH_DIRECTORY}
bed_basename=`basename $bedfile`
samplename="${bed_basename%%.*}"
cp /GPFS/zhangli_lab_permanent/maitena/GRCm39.genome.fa ${SCRATCH_DIRECTORY}/genome
cp /GPFS/Magda_lab_permanent/maitenat/scs_yth/code/metaPlotR-master.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf metaPlotR-master.tar.gz
cp /GPFS/Magda_lab_permanent/maitenat/GRCm39_gencode.genePred ${SCRATCH_DIRECTORY}

# Sort bedfile
# Select just first 6 columns of bedfile
awk -v OFS='\t' '{print $1, $2, $3, $4, $12, $6}' ${bed_basename} > file.bed 
# awk -v OFS='\t' '{print $1, $2, $3, $4, $16, $10}' ${bed_basename} > file.bed 
sort -k1,1 -k2,2n file.bed > ${samplename}.sorted.bed

# Preprocess reference genome
# Split multifasta into multiple single fastas
cd genome
awk -F " " '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=ID".fa"} {print >> F}' GRCm39.genome.fa
rm GRCm39.genome.fa
mv GL456367.1.fa chrUn_GL456372v1.fa
mv GL456381.1.fa chrUn_GL456381v1.fa
mv GL456385.1.fa chrUn_GL456385v1.fa
mv JH584304.1.fa chrUn_JH584304v1.fa
mv GL456210.1.fa chr1_GL456210v1_random.fa
cd ..
# Annotate reference genome
perl metaPlotR-master/make_annot_bed.pl --genomeDir genome/ --genePred GRCm39_gencode.genePred > GRCm39_annot.bed
sort -k1,1 -k2,2n GRCm39_annot.bed > GRCm39_sorted.annot.bed

# Catalogue transcriptome
perl metaPlotR-master/size_of_cds_utrs.pl --annot GRCm39_sorted.annot.bed > region_sizes.txt

# Annotate CIMS bedfile
perl metaPlotR-master/annotate_bed_file.pl --bed ${samplename}.sorted.bed --bed2 GRCm39_sorted.annot.bed > annot_${samplename}.bed

# Metagene coordinates
perl metaPlotR-master/rel_and_abs_dist_calc_maitena.pl --bed annot_${samplename}.bed --regions region_sizes.txt > ${samplename}.dist.measures.txt

# After the job is done we copy our output back to our directory
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}.dist.measures.txt ${outdir}
cp -r ${SCRATCH_DIRECTORY}/annot_${samplename}.bed ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
