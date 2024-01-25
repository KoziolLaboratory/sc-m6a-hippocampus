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

The `scDART` folder contains the data for the analysis of single-cell dartseq data. The starting data for these scripts are the RSEC metrics obtained from the Rhapsody software. The scripts in this section should be run in the order listed. 


1. `make_bullseye_mat.sbatch`
2. `find_edit_sites.sbatch`
3. `filter_bullseye.sbatch`
4. `integr_scDART.Rmd`
5. `scibet_tm.R`
6. `scMCA_mice.R`
7. `call_metaplotR.sbatch`
8. `m6a_analysis.Rmd`
9. `m6a_UMAP.Rmd`
10. `diff_expression_yth_mut.Rmd`

For the homogeneous m6A analysis:

1. `get_homogeneous_m6a.Rmd`
2. `call_metaplotR_homogeneous.sbatch`
3. `prepare_mk1_metaplootr.R`
4. `compare_m6a.Rmd`
5. `homogeneous_m6a_enrichment.Rmd`


### MeRIP

The `MeRIP` folder contains scripts and data related to m6A RNA immunoprecipitation (MeRIP). This technique is used to enrich m6A-modified RNAs and identify m6A sites at a transcriptome-wide level. 

For the analysis of this data, the following scripts are called, in this order:

1. `rm_rRNA.sbatch`
2. `run_cutadapt.sbatch`
3. `run_fastqc.sbatch`
4. `run_multiqc.sbatch`
5. `run_HISAT2.sbatch`
6. `peakcalling.sbatch`
7. `intersect_peaks.sbatch`
8. `annotate_genes_to_peaks.sh`
9. `venn_diagrams.Rmd`

### sc_bulk_MeRIP

The `sc_bulk_MeRIP` folder contains a script for the comparison of the results obtained by MeRIP, single-cell dartseq and bulk dartseq data. The script is `sc_bulk_MeRIP.Rmd`.



## Requirements

This project requires the following tools and libraries:

### Software
- [R](https://www.r-project.org/)
- [Perl](https://www.perl.org/)
- [Python](https://www.python.org)
- [HISAT2](http://daehwankimlab.github.io/hisat2/) for read alignment
- [BWA](http://bio-bwa.sourceforge.net/) for read alignment
- [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) for read alignment
- [samtools](http://www.htslib.org/) for manipulating aligned reads
- [HTSlib](http://www.htslib.org/) for high-throughput sequencing data processing
- [Tabix](http://www.htslib.org/doc/tabix.html) for indexing and retrieving data in tab-delimited text files
- [bedtools](https://bedtools.readthedocs.io/) for genome arithmetic
- [Picard](https://broadinstitute.github.io/picard/) for manipulating high-throughput sequencing data
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control checks on raw sequence data
- [MultiQC](https://multiqc.info/) for aggregating quality control checks
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/) for adapter trimming
- [MACS2](https://github.com/macs3-project/MACS) for peak calling
- [CTK](https://zhanglab.ccmb.med.umich.edu/CTK/) for bulk dart-seq data analysis
- [Bullseye](https://github.com/Boyle-Lab/Bullseye) for single-cell dart-seq data analysis


### R Libraries
- [tidyverse](https://www.tidyverse.org/) for data manipulation and visualization
- [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html) for generating Venn diagrams
- [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html) for color palettes
- [Seurat](https://satijalab.org/seurat/) for single cell genomics
- [Nebulosa](https://github.com/krassowski/nebulosa) for visualization of single-cell RNA-seq data
- [scibet](https://github.com/BGI-shenzhen/scibet) for cell type identification
- [scMCA](https://github.com/SCA-IRCM/SingleCellMultiOmics) for multi-omics data analysis
- [scales](https://scales.r-lib.org/) for graphical scales mapping
- [ggVolcano](https://cran.r-project.org/web/packages/ggVolcano/index.html) for creating volcano plots
- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) for comparing and visualizing functional profiles
- [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) for annotation database interfaces
- [org.Mm.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) for organism annotations
- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) for BioMart databases
- [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) for visualization of functional enrichment result
- [DOSE](https://bioconductor.org/packages/release/bioc/html/DOSE.html) for disease ontology semantic and enrichment analysis


### Data
- Reference genomes (GRCm39/mm39, hg19)
- Adapter sequences for the sequencing platform used
- Tabula Muris Brain Non-Myeloid model (GSE109774_scibet_core.csv)
- Gene annotations (GRCm39_gencode, hg19_gencode)
- Variant annotations (mgp_REL2021_snps.bed, common_all_20180423_snps.bed, GCA_000001635.9_current_ids.bed)


### Job Scheduler
- [SLURM](https://slurm.schedmd.com/overview.html) for job scheduling on a high-performance computing cluster. Many of the scripts in this project are written as SLURM batch scripts and require SLURM for execution. However, these SLURM scripts are essentially shell scripts. If you are not using SLURM, you can adapt them to run as regular bash scripts with some minimal knowledge of bash scripting.

Please ensure all dependencies are installed and the raw data is available in the appropriate directories before running the scripts.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

