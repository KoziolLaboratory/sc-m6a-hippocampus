# Single cell discovery of m6A RNA modifications in the hippocampus

This repository contains the code associated with the article titled "Single cell discovery of m6A modifications in the hippocampus." The code provided here was used for the analysis and generation of results presented in the article.

## Introduction

The code in this repository is organized into four main folders:

1. `bulkDART`: This folder contains the code used for the analysis of bulk dartseq data.

2. `scDART`: This folder contains the code used for the analysis of single-cell dartseq data.

3. `MeRIP`: This folder contains scripts related to m6A RNA immunoprecipitation (MeRIP).

4. `sc_bulk_MeRIP`: This folder contains scripts related to the comparison of the results obtained by MeRIP, single-cell dartseq and bulk dartseq results.


## Code Overview and usage

The scripts used in this work are written in Shell (`.sh`, `.sbatch`), R (`.R`, `.Rmd`), and possibly others depending on the specific requirements of the analysis. The `.sbatch` files are scripts written for SLURM job scheduler for Linux. These scripts are used to submit jobs to the SLURM scheduler, which then manages the job's execution on the cluster. 

### bulkDART

The `bulkDART` folder contains scripts for the analysis of bulk dartseq data. Some of these scripts are used in the analysis of data from the HEK293T cell line, while others are used for the analysis of data from the mice hippocampus. 

For the analysis of mice hippocampus data, the following scripts are called, in this order:

1. `run_fastqc.sbatch`
2. `rename_fq.sbatch`
3. `align_mm39.sbatch`
4. `call_CU.sbatch`
5. `filter_CIMS_mm39.sbatch`
6. `dartseq_replicate_counts_mm39.Rmd`
7. `obtain_CU_2out3.R`
8. `filter_yth.sbatch`
9. `call_metaPlotR_mm39.sbatch`
10. `hippocampus_variant_analysis.Rmd`

For the analysis of HEK293T cell line data, the following scripts are called, in this order:

1. `run_fastqc.sbatch`
2. `rename_fq.sbatch` 
3. `align_hg19.sbatch`
4. `call_CU.sbatch`
5. `filter_CIMS_hg19.sh`
6. `dartseq_replicate_counts_HEK.Rmd`
7. `obtain_CU_2out3.R`
8. `filter_yth.sbatch`
9. `call_metaPlotR_hg19.sbatch`
10. `HEK_YTH-E_vs_E-YTH.Rmd`

### scDART

TODO

### MeRIP

The `MeRIP` folder contains scripts and data related to m6A RNA immunoprecipitation (MeRIP). This technique is used to enrich m6A-modified RNAs and identify m6A sites at a transcriptome-wide level. 

(...)
n. `annotate_genes_to_peaks.sh`
n+1. `venn_diagrams.Rmd`

### sc_bulk_MeRIP

The `sc_bulk_MeRIP` folder contains a script for the comparison of the results obtained by MeRIP, single-cell dartseq and bulk dartseq data. The script is `sc_bulk_MeRIP.Rmd`.



## Requirements

...

## License

This project is licensed under the MIT License - see the LICENSE file for details.

