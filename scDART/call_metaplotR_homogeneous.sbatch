#!/bin/bash
#SBATCH --job-name=metaplotr
#SBATCH --output=metaplotr_%A.out
#SBATCH --error=metaplotr_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=80G

# Description:  Here we are calling metaPlotR on homogeneuos sites

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
cp /home/Magda_lab/maitenat/code/metaPlotR-master.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf metaPlotR-master.tar.gz
cp /GPFS/Magda_lab_permanent/maitenat/GRCm39_gencode.genePred ${SCRATCH_DIRECTORY}

# Sort bedfile
# We don't have the strand info in the file, and we need it for metaplotR (does not work putting a .), so we use the mk=1 file to get it
cell_type=${bed_basename%_mk1_common_mk09_5cells.txt}
case $cell_type in
    "immune_cells")
        strand_info_file="/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_IMC.bed"
        ;;
    "astrocytes")
        strand_info_file="/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_ASC.bed"
        ;;
    "endothelial")
        strand_info_file="/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_EC.bed"
        ;;
    "neurons")
        strand_info_file="/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_NEUR.bed"
        ;;
    "oligodendrocytes")
        strand_info_file="/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_OLG.bed"
        ;;
    *)
        echo "Error: cell type not recognized"
        ;;
esac

# print the whole lines from strand_info_file if the first three columns match the first three columns of bed_basename
awk 'NR==FNR{a[$1,$2,$3]=$0;next} ($1,$2,$3) in a{print a[$1,$2,$3]}' ${strand_info_file} ${bed_basename} > file.bed
sort -k1,1 -k2,2n file.bed > ${samplename}.sorted.bed

# break if the number of lines in ${samplename}.sorted.bed is not 1 less than the number of lines in ${bed_basename}
# minus1 is because bedfile has a header
num_lines_bed_basename=$(wc -l < "${bed_basename}")
num_lines_sorted_bed=$(wc -l < "${samplename}.sorted.bed")

if [ $((num_lines_bed_basename - 1)) -ne $num_lines_sorted_bed ]; then
    echo "Error: Number of lines in ${samplename}.sorted.bed is not 1 less than the number of lines in ${bed_basename}"
    exit 1
fi

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
