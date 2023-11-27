#!/bin/bash
#SBATCH --job-name=bullseye
#SBATCH --output=bullseye_edsite_%A.out
#SBATCH --error=bullseye_edsite_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=80G
#SBATCH --partition=q_fat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maitenat@gmail.com

module load bullseye/bullseye_20220908

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

# Copy necessary files to scratch directory
mat_filename=$1
refflat=refFlat_ncbi-refseq_mm39.txt
REF=GRCm39.genome.fa
cp $mat_filename ${SCRATCH_DIRECTORY}
cp $mat_filename".tbi" ${SCRATCH_DIRECTORY}
cp $refflat ${SCRATCH_DIRECTORY}
cp $REF ${SCRATCH_DIRECTORY}
refflat_basename=`basename $refflat`
mat_basename=`basename ${mat_filename}`
samplename="${mat_basename%%.*}"
ref_basename=`basename ${REF}`

echo 'processing' $samplename

# Run find_edit_site.pl to identify m6A sites
echo 'identifying m6a sites'

perl /usr/nzx-cluster/apps/bullseye/Bullseye/Code/Find_edit_site.pl --annotationFile $refflat_basename --EditedMatrix ${mat_filename} --editType C2U --barcode --genome $ref_basename --minEdit 5 --maxEdit 100 --editFoldThreshold 1.5 --MinEditSites 2 --EditedMinCoverage 10 --ControlMinCoverage 10 --cpu 8 --outfile ${samplename}_edit-sites_bc.bed --verbose

# After the job is done we copy our output back to our directory
outdir=/home/Magda_lab/maitenat/scratch60/scs_yth/${samplename}/bullseye/
mkdir -p ${outdir} && cp -r ${samplename}_edit-sites_bc.bed ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
