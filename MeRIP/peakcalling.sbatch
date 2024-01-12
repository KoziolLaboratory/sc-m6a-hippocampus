#!/bin/bash
#SBATCH -J macs2_peakcalling
#SBATCH -o macs2_peakcalling.out
#SBATCH -p q_cn
#SBATCH -n 1

module load macs2/2.2.6

# Parameters of the current job
pwd; hostname; date

bamfile=$1
samplename=`echo $bamfile | sed 's/\.Hisat_aln\.sorted\.bam//'`
control_bamfile=$2
outdir=/home/Magda_lab/maitenat/scratch60/MeRIP/aligned/peakcalling

# Copy necessary files to scratch directory
cp $bamfile ${SCRATCH_DIRECTORY}
cp $control_bamfile ${SCRATCH_DIRECTORY}
cp ${bamfile}.bai ${SCRATCH_DIRECTORY}
cp ${control_bamfile}.bai ${SCRATCH_DIRECTORY}
bamfile_basename=`basename $bamfile`
control_bamfile_basename=`basename $control_bamfile`

macs2 callpeak -t $bamfile_basename -c $control_bamfile_basename -B --SPMR --keep-dup all -n $samplename --outdir $outdir -f BAM -g 105516336 --nomodel --extsize 150 -q 0.05 --fe-cutoff 2


