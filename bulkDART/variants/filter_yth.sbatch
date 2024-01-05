#!/bin/bash
#SBATCH --job-name=filt_yth
#SBATCH --output=filter_yth_%A.out
#SBATCH --error=filter_yth_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=60G

# Description:  Here we filter the CU variants in yth and ythmut samples. Specifically we remove APOBEC and WT 
# from YTH and YTH-MUT. From YTH, we keep only those variants with edition frequency 1.5X higher than YTH-MUT sample

module load bedtools/2.29.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
yth=$1 
ythmut=$2 
wt=$3 
apobec=$4 
samplename=$5

# Copy necessary files to scratch directory
cp $yth ${SCRATCH_DIRECTORY}
cp $ythmut ${SCRATCH_DIRECTORY}
cp $wt ${SCRATCH_DIRECTORY}
cp $apobec ${SCRATCH_DIRECTORY}

yth_basename=`basename ${yth}`
ythmut_basename=`basename ${ythmut}`
wt_basename=`basename ${wt}`
apobec_basename=`basename ${apobec}`

# Discard wt and apobec
## yth
bedtools subtract -s -A -a ${yth_basename} -b ${wt_basename} > yth_discard1.cims
bedtools subtract -s -A -a yth_discard1.cims -b ${apobec_basename} > yth_discard.cims
## yth-mut
bedtools subtract -s -A -a ${ythmut_basename} -b ${wt_basename} > ythmut_discard1.cims
bedtools subtract -s -A -a ythmut_discard1.cims -b ${apobec_basename} > ythmut_discard.cims

# For yth, keep only sites with 1.5X enrichment over yth-mut
### For that, first, find C to T transitions that are present only in APOBEC1-YTH
bedtools subtract -s -A -a yth_discard.cims -b ythmut_discard.cims > yth_notythmut.cims
## Next, merge the YTH and mut files for common sites.
bedtools intersect -s -wa -wb -a yth_discard.cims -b ythmut_discard.cims > yth_mergeythmut.cims
## For those, identify sites that have a m/k ratio that is 1.5-fold greater than in APOBEC1-YTHmut samples:
awk '{a=$7/$14;print $0,a;}' yth_mergeythmut.cims | awk '$15 >= 1.5 {print $0}' > yth_mk1.5overythmut.cims
## Add back the sites in YTH that were not in mut:
cat yth_mk1.5overythmut.cims yth_notythmut.cims > ${samplename}_yth_filt-all.cims

# After the job is done we copy our output back to our directory
outdir=/GPFS/Magda_lab_temp/maitenat/dartseq_replicates/aligned/
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}_yth_filt-all.cims ${outdir}
mv yth_discard.cims ${samplename}_yth_no-back-apo.cims
mv ythmut_discard.cims ${samplename}_ythmut_no-back-apo.cims
cp -r ${SCRATCH_DIRECTORY}/${samplename}_yth_no-back-apo.cims ${outdir}
cp -r ${SCRATCH_DIRECTORY}/${samplename}_ythmut_no-back-apo.cims ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
