#!/bin/bash
#SBATCH --job-name=bullseye
#SBATCH --output=bullseye_%A.out
#SBATCH --error=bullseye_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=80G
#SBATCH --partition=q_fat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maitenat@gmail.com


module load bullseye/bullseye_20220908
module load htslib/1.9
module load tabix/0.2.6

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

# Copy necessary files to scratch directory
bam_filename=$1
barcode_filename=$2

refflat=refFlat_ncbi-refseq_mm39.txt
REF=GRCm39.genome.fa
cp $bam_filename ${SCRATCH_DIRECTORY}
cp $bam_filename".bai" ${SCRATCH_DIRECTORY}
cp $refflat ${SCRATCH_DIRECTORY}
bam_basename=`basename $bam_filename`
refflat_basename=`basename $refflat`

# we need to remove uppercase from the bamfile extension for bullseye not to break
samplename="${bam_basename%%.*}"
short_samplename="${bam_basename%%_final*}"
mv $bam_basename $samplename".bam"
# same with the index file
mv $bam_basename".bai" $samplename".bam.bai"

echo 'processing' $samplename".bam"

# Run parseBAM.pl to count the nucleotides at each position in the genome
perl /usr/nzx-cluster/apps/bullseye/Bullseye/Code/parseBAM.pl --mode SingleCell --Cell_ID_pattern CB:Z: --input $samplename".bam" --cpu 10  --minCoverage 10 --filterBarcodes $barcode_filename --output ${short_samplename}.matrix --verbose


# After the job is done we copy our output back to our directory
outdir=/home/Magda_lab/maitenat/scratch60/scs_yth/${short_samplename}/bullseye/
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${short_samplename}.matrix.gz ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}