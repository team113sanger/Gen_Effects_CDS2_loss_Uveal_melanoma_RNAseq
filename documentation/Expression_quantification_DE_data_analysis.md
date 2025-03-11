# Differential expression analysis

## Overview

This document provides the code and rationale of the approach used to identify differentially expressed genes between the multiple paired comparisons made amongst the different experimental groups on both _CDS2_ targeted SW837 cell-lines or OMM2.5 xenografts.

## Required dependencies

### Required software

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.2.2`[**here**](https://cran.r-project.org/)
- **sleuth**: sleuth `0.30.1` [**here**](https://github.com/pachterlab/sleuth/tree/v0.30.1) - Installed with renv
- **fgsea**: sleuth `1.30.0` [**here**](https://bioconductor.org/packages/release/bioc/html/fgsea.html) - Installed with renv
- **IGV**: IGV `2.13.1` [**here**](https://data.broadinstitute.org/igv/projects/downloads/2.13/)

### R environment setup for sleuth
To reproduce the R environment used to generate the figures install the packages required to generate the final files and plots, the following steps are required run the following commands within `R v4.4.1`. 

Set the `projectdir` variable to the path where the repository was cloned into and run the following commands in R to install the required packages and restore the environment from the `renv.lock` file.

```R
# Run in R v4.4.1
projectdir<-"/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES"
plotsc_dir<- file.path(projectdir,"scripts/expression")

setwd(plotsc_dir)
install.packages("renv")
library(renv)
# To rebuild an environment from the renv.lockfile
# use renv::restore(lockfile=file.path(plotsc_dir, "renv.lock"), prompt=FALSE) if you want R not to prompt you for confirmation
renv::restore(lockfile=file.path(plotsc_dir, "renv.lock"))

```

## Result reproduction

### Differential expression analyses

To identify differential expresses genes in _CDS2_ targeted lines a total of 5 paired comparisons were made based on the `Groups` column on the **7687_3364_sample_metadata.tsv** table. 
**Two** for the experiments with SW837 cells:
- **SW837_C9_safe_gRNA VS SW837_CDS2**: Comparison between **_CDS2_ targeted SW837** lines and safe targeting gRNA (_STG1_) SW837-Cas9 targeted cells.
- **SW837_uninfected VS SW837_CDS2**: Comparison between **_CDS2_ targeted SW837** lines and uninfected SW837-Cas9 cells.

**Three** for the OMM2.5 grafts:
- **OMMCDS2_Tum_Cas9WT VS CDS2_Tum_Treated**: Comparison between **_CDS2_ targeted OMM2.5** treated with doxocycline and OMM2.5-Cas9 treated with doxocycline but no gRNA.
- **OMMCDS2_Tum_UnTreated VS CDS2_Tum_Treated**: Comparison between **_CDS2_ targeted OMM2.5** without doxocycline and OMM2.5-Cas9 treated with inducible _CDS2_ gRNA construct without doxocycline.
- **OMMCDS2_Tum_Cas9WT VS CDS2_Tum_UnTreated**: Comparison between OMM2.5-Cas9 grafts treated with doxocycline but no gRNA and **_CDS2_ targeted OMM2.5** without doxocycline.

Three samples were removed from OMM2.5 analysis due to higher Jensen-Shannon divergence  and lower Pearson correlation than the rest of the samples within their group. The samples were: CDS2_Tumour_Cas9_WT_R4_RNA, CDS2_Tumour_R2_RNA and CDS2_Tumour_Treated_R1_RNA

Significant differentially expressed genes sets identified with qval<0.01 are then used to run a gene set enrichment analysis using GSEA to assess the gene set with three of the [Human Molecular Signatures Database (MSigDB) gene sets](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp): H:Hallmark, C2:CP - Curated gene sets Canonical Pathways and C6: oncogenic signatures. 


#### Required input files 

To be able to run the the analysis with sleuth, the following files are required to be present within the repository:

- **`kallisto`**: This is the folder with all the results per sample of the `kallisto quant`. This can be downloaded from **figshare** [**here**](https://figshare.com/ndownloader/articles/28016804/versions/1?folder_path=kallisto)
- [**7687_3364_sample_metadata.tsv**](../metadata/7687_3364_sample_metadata.tsv) : This is the metadata file containing the information of the samples used in the analysis, treatment, and targeting status. This file is located in the `metadata` directory of the repository.

#### Execution of `DE_data_analysis_and_TPM_table_collation.R `

The analysis below was run in a machine with 7 threads and 16GB of RAM. The script has some values by default but can be modified if a shell environment variable is present, see below:
 - **DE_SEED**: The value of the seed used to ensure reproducibility was **42**, this can be changed by setting the shell environment variable **DE_SEED**
- **DE_NCORES**: The number of threads used by default is 7, this can be changed by setting the shell environment variable **DE_NCORES** to something different. 
- **DE_RESULTS_DIR**: The default location of the results directory is `analysis/expression_analysis`, however this can be set by providing setting the shell environment variable **DE_RESULTS_DIR**

```bash
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq
STUDY=7687

# Get the kallisto results from figshare
cd ${PROJECTDIR}/analysis/
wget -c "https://figshare.com/ndownloader/articles/28016804/versions/1?folde
r_path=kallisto" -O kallisto.zip
unizp ${PROJECTDIR}/analysis/kallisto.zip
rm ${PROJECTDIR}/analysis/kallisto.zip


cd ${PROJECTDIR}/scripts/expression/
Rscript ${PROJECTDIR:?unset}/scripts/expression/DE_data_analysis_and_TPM_table_collation.R 


```

## Results

### IGV images 

The IGV images showing _CDS2_ targeting sites used in the paper are located in the [`analysis/IGV_RNAseq`](../analysis/IGV_RNASeq/) directory. The images are named according to the experiment cell-line displayed, **SW837** or **OMM2.5**. Only Human BAM files and filtered Human BAM files were respectively used.

Files used in paper figures:
- [Extended Data Fig18](../analysis/IGV_RNASeq/CDS2_Xfilt_RNASEq_ALL_samples_CDS2_ontar_site_10bp_15bp_expanded.png)


### Expression analysis results

The list of files, tables and plots generated by the code above are created in the `analysis/expression_analysis` directory by default. However given the number of files and file sizes, the results were uploaded to **figshare** and can be downloaded from the **figshare** dataset [Expression data and differential expression analysis of CDS2 targeted SW837 and OMM2.5 xenografted cells](https://figshare.com/articles/dataset/Expression_data_and_differential_expression_analysis_of_i_CDS2_i_targeted_SW837_and_OMM2_5_xenografted_cells/28016804). The files are located within the folder with the same name `expression_analysis`.

A total of three different directories are created:
- **cohort_tables**: This folder contains the tables with the matrices containing the collated counts or Transcripts per Million (TPM) values for all the genes int he annotation for all samples.

- **DE_results**: This folder contains a folder, one for each of the five different comparisons. Each folder contains the results of the differential expression analysis for the respective paired comparison, this includes tables with set of genes that were found to be differentially expressed the matrices containing the collated counts or Transcripts per Million (TPM) values for all the genes int he annotation for all samples. Every comparison folder has a compressed `.RDs` file which contains the `sleuth_object` with the results of the comparison and parameters which can be reloaded into R. The list of files generated have the name of the sample comparison as the `<PREFIX>` of the file names e.g. `SW837_C9_safe_gRNA_VS_SW837_CDS2`. The description of the files generate for each comparison are:
    - **Lrt_gsea** : Folder with the GSEA results for the DE gene sets found using sleuth's likelihood ratio test
	- **`<PREFIX>`_C2_CP_GSEA_padj_0.01box.pdf**: box plot with C2:CP GSEA results for the DE gene sets found using sleuth wald test
	- **`<PREFIX>`_C2_CP_GSEA_padj_0.01.tsv**: table with C2:CP GSEA results for DE genes found using wald test
	- **`<PREFIX>`_C6_Cancer_sig_GSEA_padj_0.01box.pdf**: box plot with C6 GSEA results for the DE gene sets found 
	- **`<PREFIX>`_C6_Cancer_sig_GSEA_padj_0.01.tsv**: table with C2:CP GSEA results for DE genes found using wald test
	- **`<PREFIX>`_CDS2_CHAC1_DDIT4_SESN2_TPM_pergrp.pdf** : box plots with the observed TPM expression in _CDS2,CHAC1,DDIT4 and SESN2_ genes  
	- **`<PREFIX>`_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf**: boxplot showing the expression of _CDS2_ based on scaled reads
	- **`<PREFIX>`_CDS2_ENSG00000101290_TPM_pergrp.pdf**: boxplot showing the expression of _CDS2_ based on TPM
	- **`<PREFIX>`_DE_LRT_sig_0.01_perGeneRes.tsv**: table with the list of DE genes found using the likelihood ratio test
	- **`<PREFIX>`_DE_LRT_sig_0.01.tsv**: table with the list of DE transcripts found using the likelihood ratio test
	- **`<PREFIX>`_DE_WaldTest_sig_0.01_perGeneRes_log2FCFilt.tsv**: table with the list of DE genes found using the wald test and absolut(log2FC)>0.5849625
	- **`<PREFIX>`_DE_WaldTest_sig_0.01_perGeneRes.tsv**: table with the list of DE genes found using the wald test and absolut(log2FC)>0.5849625
	- **`<PREFIX>`_FRRS1L_ENSG00000260230_TPM_pergrp.pdf**: boxplot showing the expression of _FRRS1L_ based on TPM a potential off-target from _CDS2_ gRNA used
	- **`<PREFIX>`_Hallmark_GSEA_padj_0.01box.pdf**: boxplot Hallmark GSEA results for DE genes found using wald test 
	- **`<PREFIX>`_Hallmark_GSEA_padj_0.01.tsv**: table with Hallmark GSEA results for DE genes found using wald test
	- **`<PREFIX>`_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf**: heatmap with the Jensen Shannon divergence
	- **`<PREFIX>`_MA_plot_WaldTest_treatment.pdf**: MA plot showing DE genes found using the wald test
	- **`<PREFIX>`_PCA_scaledreads_per_base_varexp.pdf**: PCA showing PC1 and PC2 based on the scaled reads per gene, as well as barplot with percentage of variance explained per PC.
	- **`<PREFIX>`_PLD4_ENSG00000166428_TPM_pergrp.pdf**: boxplot showing the expression of _PLD4_ based on TPM
	- **`<PREFIX>`_QQ_plotPCA_scaledreads_per_base_lrt.pdf**: QQplot based on scaled read expression values
	- **`<PREFIX>`_scaledreads_per_base_distribution.pdf**: histogram showing the distribution of expression values observed based on log(mean(scaled_reads_per_base)+1)
	- **`<PREFIX>`_sleuth_object.RDS**: Compressed `.Rds` file that contains the Sleuth object with all the results and analysis for this comparison.
	- **`<PREFIX>`_SMG1_ENSG00000157106_TPM_pergrp.pdf:** boxplot showing the expression of _SMG1_ based on TPM a potential off-target from _CDS2_ gRNA used
	- **`<PREFIX>`_TPM_distribution.pdf:** histogram showing the distribution of expression values per group based on log(mean(TPM)+1)
	- **`<PREFIX>`_Volcano_treatment.pdf:** Volcano plot showing bet values(x axis) vs -log10(qval) obtained for DE genes found using the Wald test
    - **`<PREFIX>`_Volcano_treatment.png:** Volcano plot showing bet values(x axis) vs -log10(qval) obtained for DE genes found using the Wald test
- **overall_qc_figures**: plots showing different heatmap plots with the Jensen-Shannon divergence and Pearson correlation across all the possible paired comparisons of the samples, calculated from the log2(TPM+1) values of expression.

**Files used in the paper figures:**
- [`SW837_C9_safe_gRNA_VS_SW837_CDS2_CDS2_ENSG00000101290_TPM_pergrp.pdf`](https://figshare.com/ndownloader/files/52504022): was used in **Extended Data Fig12**
- [`SW837_C9_safe_gRNA_VS_SW837_CDS2_DE_WaldTest_sig_0.01_perGeneRes_log2FCFilt.tsv`](https://figshare.com/ndownloader/files/52504034): was used in **Extended Data Fig12** legend

#### Summary of files contained 

```bash
PROJECTDIR=
cd ${PROJECTDIR}/analysis/expression_analysis
tree -L 3
├── cohort_tables
│   ├── 7687_kallisto_est_count_ENSv103.tsv.gz
│   └── 7687_kallisto_est_tpm_ENSv103.tsv.gz
├── DE_results
│   ├── CDS2_targeting_kallisto_sleuth_env.RData
│   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated
│   │   ├── Lrt_gsea
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_C2_CP_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_C2_CP_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_C6_Cancer_sig_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_C6_Cancer_sig_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_CDS1_ENSG00000163624_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_CDS2_ENSG00000101290_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_CRX_ENSG00000105392_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_DE_LRT_sig_0.01_perGeneRes.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_DE_LRT_sig_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_DE_WaldTest_sig_0.01_perGeneRes_log2FCFilt.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_DE_WaldTest_sig_0.01_perGeneRes.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_FRRS1L_ENSG00000260230_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_Hallmark_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_Hallmark_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_MA_plot_WaldTest_treatment.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_NRL_ENSG00000129535_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_OSGIN1_ENSG00000140961_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_PCA_scaledreads_per_base_varexp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_QQ_plotPCA_scaledreads_per_base_lrt.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_scaledreads_per_base_distribution.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_sleuth_object.RDS
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_SMG1_ENSG00000157106_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_SNX19_ENSG00000120451_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_TPM_distribution.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_Volcano_treatment.pdf
│   │   └── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated_Volcano_treatment.png
│   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated
│   │   ├── Lrt_gsea
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_C2_CP_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_C2_CP_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_C6_Cancer_sig_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_C6_Cancer_sig_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_CDS1_ENSG00000163624_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_CDS2_ENSG00000101290_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_DE_LRT_sig_0.01_perGeneRes.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_DE_LRT_sig_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_DE_WaldTest_sig_0.01_perGeneRes_log2FCFilt.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_DE_WaldTest_sig_0.01_perGeneRes.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_FRRS1L_ENSG00000260230_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_Hallmark_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_Hallmark_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_MA_plot_WaldTest_treatment.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_OSGIN1_ENSG00000140961_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_PCA_scaledreads_per_base_varexp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_QQ_plotPCA_scaledreads_per_base_lrt.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_scaledreads_per_base_distribution.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_sleuth_object.RDS
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_SMG1_ENSG00000157106_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_SNX19_ENSG00000120451_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_TPM_distribution.pdf
│   │   ├── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_Volcano_treatment.pdf
│   │   └── OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated_Volcano_treatment.png
│   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated
│   │   ├── Lrt_gsea
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_C2_CP_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_C2_CP_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_C6_Cancer_sig_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_C6_Cancer_sig_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_CDS1_ENSG00000163624_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_CDS2_ENSG00000101290_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_DE_LRT_sig_0.01_perGeneRes.tsv
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_DE_LRT_sig_0.01.tsv
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_DE_WaldTest_sig_0.01_perGeneRes_log2FCFilt.tsv
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_DE_WaldTest_sig_0.01_perGeneRes.tsv
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_FRRS1L_ENSG00000260230_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_Hallmark_GSEA_padj_0.01box.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_Hallmark_GSEA_padj_0.01.tsv
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_MA_plot_WaldTest_treatment.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_OSGIN1_ENSG00000140961_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_PCA_scaledreads_per_base_varexp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_QQ_plotPCA_scaledreads_per_base_lrt.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_scaledreads_per_base_distribution.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_sleuth_object.RDS
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_SMG1_ENSG00000157106_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_SNX19_ENSG00000120451_TPM_pergrp.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_TPM_distribution.pdf
│   │   ├── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_Volcano_treatment.pdf
│   │   └── OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated_Volcano_treatment.png
│   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2
│   │   ├── Lrt_gsea
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_C2_CP_GSEA_padj_0.01box.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_C2_CP_GSEA_padj_0.01.tsv
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_C6_Cancer_sig_GSEA_padj_0.01box.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_C6_Cancer_sig_GSEA_padj_0.01.tsv
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_CDS2_CHAC1_DDIT4_SESN2_TPM_pergrp.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_CDS2_ENSG00000101290_TPM_pergrp.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_DE_LRT_sig_0.01_perGeneRes.tsv
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_DE_LRT_sig_0.01.tsv
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_DE_WaldTest_sig_0.01_perGeneRes_log2FCFilt.tsv
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_DE_WaldTest_sig_0.01_perGeneRes.tsv
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_FRRS1L_ENSG00000260230_TPM_pergrp.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_Hallmark_GSEA_padj_0.01box.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_Hallmark_GSEA_padj_0.01.tsv
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_MA_plot_WaldTest_treatment.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_PCA_scaledreads_per_base_varexp.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_PLD4_ENSG00000166428_TPM_pergrp.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_QQ_plotPCA_scaledreads_per_base_lrt.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_scaledreads_per_base_distribution.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_sleuth_object.RDS
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_SMG1_ENSG00000157106_TPM_pergrp.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_TPM_distribution.pdf
│   │   ├── SW837_C9_safe_gRNA_VS_SW837_CDS2_Volcano_treatment.pdf
│   │   └── SW837_C9_safe_gRNA_VS_SW837_CDS2_Volcano_treatment.png
│   └── SW837_uninfected_VS_SW837_CDS2
│       ├── Lrt_gsea
│       ├── SW837_uninfected_VS_SW837_CDS2_C2_CP_GSEA_padj_0.01box.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_C2_CP_GSEA_padj_0.01.tsv
│       ├── SW837_uninfected_VS_SW837_CDS2_C6_Cancer_sig_GSEA_padj_0.01box.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_C6_Cancer_sig_GSEA_padj_0.01.tsv
│       ├── SW837_uninfected_VS_SW837_CDS2_CDS2_CHAC1_DDIT4_SESN2_TPM_pergrp.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_CDS2_ENSG00000101290_TPM_pergrp.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_DE_LRT_sig_0.01_perGeneRes.tsv
│       ├── SW837_uninfected_VS_SW837_CDS2_DE_LRT_sig_0.01.tsv
│       ├── SW837_uninfected_VS_SW837_CDS2_DE_WaldTest_sig_0.01_perGeneRes_log2FCFilt.tsv
│       ├── SW837_uninfected_VS_SW837_CDS2_DE_WaldTest_sig_0.01_perGeneRes.tsv
│       ├── SW837_uninfected_VS_SW837_CDS2_FRRS1L_ENSG00000260230_TPM_pergrp.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_Hallmark_GSEA_padj_0.01box.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_Hallmark_GSEA_padj_0.01.tsv
│       ├── SW837_uninfected_VS_SW837_CDS2_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_MA_plot_WaldTest_treatment.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_PCA_scaledreads_per_base_varexp.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_QQ_plotPCA_scaledreads_per_base_lrt.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_scaledreads_per_base_distribution.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_sleuth_object.RDS
│       ├── SW837_uninfected_VS_SW837_CDS2_SMG1_ENSG00000157106_TPM_pergrp.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_TPM_distribution.pdf
│       ├── SW837_uninfected_VS_SW837_CDS2_Volcano_treatment.pdf
│       └── SW837_uninfected_VS_SW837_CDS2_Volcano_treatment.png
└── overall_qc_figures
    ├── ALL_OMM2.5_CDS2_ENSG00000101290_TPM_pergrp.pdf
    ├── ALL_OMM2.5_CDS2_ENSG00000101290_TPM_pergrp_woFsampl.pdf
    ├── ALL_OMM2.5_CDS2reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
    ├── ALL_OMM2.5_CDS2reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base_woFsampl.pdf
    ├── ALL_Samples_reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
    ├── ALL_SW837_reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf
    ├── TPM_kallisto_gene_Pcor_heatmap_7687_samples.pdf
    ├── TPM_kallisto_gene_Pcor_heatmap_7687_samples_pGroups.pdf
    ├── TPM_kallisto_transcripts_Pcor_heatmap_7687_samples.pdf
    └── TPM_kallisto_transcripts_Pcor_heatmap_7687_samples_pGroups.pdf

13 directories, 142 files
```


