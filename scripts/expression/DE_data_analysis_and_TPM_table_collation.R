#!/usr/bin/env Rscript
############################
# DE_data_analysis_and_TPM_table_collation.R
#
#
#
#
#
############################
default_seed<-42
default_ncores<-7
envvar_seed<-"DE_SEED"
envvar_ncores<-"DE_NCORES"
envvar_results_dir<-"DE_RESULTS_DIR"


progname<-"DE_data_analysis_and_TPM_table_collation.R"

library(here)
library(readr)
library(sleuth)
library(ggplot2)
library(biomaRt)
library(RColorBrewer)
library(gplots)
# GSEA
library(fgsea)
library(msigdbr)

# Variables with the directories 
projdir<-here()
setwd(projdir)
studyid<-7687
ncores<-as.numeric(Sys.getenv(envvar_ncores, unset=default_ncores))
kalli_dir<-file.path(projdir, "analysis", "kallisto")
default_results_dir<-file.path(projdir, "analysis", "expression_analysis")
results_dir<-Sys.getenv(envvar_results_dir, unset=default_results_dir)  # Optional env-var to specify the results dir -- used for testing
tables_dir<-file.path(results_dir, "cohort_tables") 
ov_fig_dir<-file.path(results_dir, "overall_qc_figures") 
DE_res_dir<-file.path(results_dir, "DE_results")
metadata_fname<-file.path(projdir, "metadata", "7687_3364_sample_metadata.tsv")

dir.create(tables_dir, recursive = TRUE)
dir.create(ov_fig_dir, recursive = TRUE)
dir.create(DE_res_dir, recursive = TRUE)

#############################  For Sleuth USE - Get the transcript information and gene ID table ----
human<- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://feb2021.archive.ensembl.org")
#human_v112<- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://may2024.archive.ensembl.org")
tx2gene <- function(){
  hmart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://feb2021.archive.ensembl.org")
  #hmart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://may2024.archive.ensembl.org")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name","ensembl_transcript_id_version", 
                                       "gene_biotype"), mart = hmart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id_version,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}

t2g <- tx2gene()
ensv_genes<-biomaRt::getBM(attributes= c( "ensembl_gene_id", "hgnc_symbol", "hgnc_id","external_gene_name","external_gene_source",
                                            "gene_biotype", "chromosome_name","start_position","end_position"),
                              filters="ensembl_gene_id", values=unique(t2g$ens_gene),
                              mart = human)
# Get the manifest with the metadata
all(dir.exists(kalli_dir))
mdata<-readr::read_delim(metadata_fname)
samples<-mdata$RNA_Sample_name
kallisto_dirs<-sapply(samples, function(ID) file.path(kalli_dir, ID))

#Red the metadata into a table to prep the counts
s2c <- data.frame(path=kallisto_dirs, sample=samples, CDS_groups = mdata$Group, stringsAsFactors=FALSE)
s2c_cds2<- s2c[s2c$CDS_groups %in% c("CDS2_Tumour_Cas9_WT", "CDS2_Tumour_Non-treated", "CDS2_Tumour_Treated"), ]
s2c_sw<-s2c[s2c$CDS_groups %in% c("SW837_C9_CDS2_gRNA", "SW837_C9_safe_gRNA", "SW837_C9_uninfected" ), ]


# CDS2 OMM2.5 lines - Full set
suppressWarnings(so_cds2 <- sleuth_prep(s2c_cds2, ~CDS_groups, target_mapping = t2g, 
                  aggregation_column = "ens_gene", 
                  read_bootstrap_tpm=TRUE,
                  gene_mode= TRUE, 
                  extra_bootstrap_summary=TRUE, num_cores = ncores )
)
# PLoT gene Jensen-Shannon divergence per samples
print("Generating Jensen_Shannon_divergence plot ...", stdout())
pdf(file.path(ov_fig_dir, paste("ALL_OMM2.5_CDS2reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf", sep="")),
    height = 8, width = 12)
plot_sample_heatmap(obj = so_cds2)
dev.off()
pdf(file.path(ov_fig_dir, paste("ALL_OMM2.5_CDS2_ENSG00000101290_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = so_cds2, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
dev.off()

# CDS2 OMM2.5 lines - Without failed samples
failed_samples<-c("CDS2_Tumour_Cas9_WT_R4_RNA","CDS2_Tumour_R2_RNA", "CDS2_Tumour_Treated_R1_RNA")
suppressWarnings(so_cds2_wofsamp <- sleuth_prep(s2c_cds2[!(s2c_cds2$sample %in% failed_samples),  ], ~CDS_groups, target_mapping = t2g, 
                                        aggregation_column = "ens_gene", 
                                        gene_mode= TRUE, 
                                        extra_bootstrap_summary=TRUE,
                                        read_bootstrap_tpm=TRUE,
                                        num_cores = ncores )
)
# PLoT gene Jensen-Shannon divergence per samples
print("Generating Jensen_Shannon_divergence plot ...", stdout())
pdf(file.path(ov_fig_dir, paste("ALL_OMM2.5_CDS2reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base_woFsampl.pdf", sep="")),
    height = 8, width = 12)
plot_sample_heatmap(obj = so_cds2_wofsamp)
dev.off()
pdf(file.path(ov_fig_dir, paste("ALL_OMM2.5_CDS2_ENSG00000101290_TPM_pergrp_woFsampl.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = so_cds2_wofsamp, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
dev.off()

# SW837 line
suppressWarnings(
so_sw <- sleuth_prep(s2c_sw, ~CDS_groups, target_mapping = t2g, 
                      aggregation_column = "ens_gene", 
                      gene_mode= TRUE, 
                      extra_bootstrap_summary=TRUE ,num_cores = ncores)
)
# PLoT gene Jensen-Shannon divergence per samples
print("Generating Jensen_Shannon_divergence plot ...", stdout())
pdf(file.path(ov_fig_dir, paste("ALL_SW837_reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf", sep="")),
    height = 8, width = 12)
plot_sample_heatmap(obj = so_sw)
dev.off()

# Get the transcript matrix for the entire cohort with all the transcript TPM per sample ----
#READ THE SAMPLES
master_table_trans_tpm<-NULL
master_table_trans_count<-NULL
sample_list<-samples
system.time(for(i in 1:length(sample_list)){
  temp<-NULL
  #Read the sample
  temp<-data.table::fread(file=file.path(kalli_dir,sample_list[i], "abundance.tsv.gz"), header=TRUE, sep="\t", stringsAsFactors = F )
  temp<-as.data.frame(temp)
  if(i==1){
    # For Counts
    master_table_trans_count<-cbind(temp[,"est_counts"])
    rownames(master_table_trans_count)<-temp$target_id #Setting up the ENSEMBL Transcript IDs
    colnames(master_table_trans_count)<-sample_list[i]
    # For TPM
    master_table_trans_tpm<-cbind(temp[,"tpm"])
    rownames(master_table_trans_tpm)<-temp$target_id #Setting up the ENSEMBL Transcript IDs 
    colnames(master_table_trans_tpm)<-sample_list[i]
    
  }else{
    # For Counts
    master_table_trans_count<-cbind(master_table_trans_count, temp[,"est_counts"])
    colnames(master_table_trans_count)<-sample_list[1:i]
    rownames(master_table_trans_count)<-temp$target_id #Setting up the ENSEMBL Transcript IDs 
    # For TPM
    master_table_trans_tpm<-cbind(master_table_trans_tpm, temp[,"tpm"])
    colnames(master_table_trans_tpm)<-sample_list[1:i]
    rownames(master_table_trans_tpm)<-temp$target_id #Setting up the ENSEMBL Transcript IDs 
  }
}
)
# Write the transcript tables 
ensembl_transcript_id<-master_table_trans_count[,1]
data.table::fwrite(cbind(ensembl_transcript_id,master_table_trans_count[,2:dim(master_table_trans_count)[2]]), 
                   file = file.path(tables_dir, paste0(studyid, "_kallisto_est_count_ENSv112.tsv.gz" )), quote=FALSE, row.names = FALSE, sep="\t" )
ensembl_transcript_id<-master_table_trans_tpm[,1]
data.table::fwrite(cbind(ensembl_transcript_id,master_table_trans_tpm[,2:dim(master_table_trans_count)[2]]), 
                   file = file.path(tables_dir, paste0(studyid, "_kallisto_est_tpm_ENSv112.tsv.gz" )), quote=FALSE, row.names = FALSE, sep="\t" )


# Use Sleuth to read in All the samples 
suppressWarnings(so_all<-sleuth_prep(s2c, ~CDS_groups, target_mapping = t2g, 
                    aggregation_column = "ens_gene", 
                    #                       read_bootstrap_tpm=TRUE, 
                    gene_mode=TRUE,
                    extra_bootstrap_summary=TRUE, num_cores=ncores
                  ))

# PLoT gene Jensen-Shannon divergence per samples
print("Generating Jensen_Shannon_divergence plot ...", stdout())
pdf(file.path(ov_fig_dir, paste("ALL_Samples_reps_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf", sep="")),
    height = 8, width = 12)
plot_sample_heatmap(obj = so_all)
dev.off()

# Transform the tidy like table to a matrix with ENS_gene_IDs on the rows and  Samples per Column
tpm_tab<-reshape2::dcast(so_all$obs_norm, target_id ~ sample, value.var = "tpm")
dim(tpm_tab)

#  Transcripts humane the HEATMAP for the Pearson correlation comparison of all sample based on TPMs  ------------------
x<-as.matrix(master_table_trans_tpm[,2:dim(master_table_trans_count)[2]] )
x<-as.matrix(x[,2:dim(x)[2]])
#Calculate the LOG2
x<-log2(x+1)
#Calculate the pearson Correlation coefficient
xc<-cor(x)
hmcolors<-colorRampPalette(rev(brewer.pal(n=11,name = "RdYlBu" )))(60)

#GET the heatmap
pdf(file.path(ov_fig_dir, paste("TPM_kallisto_transcripts_Pcor_heatmap_", studyid,"_samples.pdf", sep="")),
    height = 12, width = 12)
gplots::heatmap.2(xc, col=hmcolors, main="Pearson correlation based on all\n transcripts \n log2(TPM+1)", 
          key=TRUE, keysize=1.5,  cexCol=0.5, cexRow=0.5, trace="none", margins=c(7,7))
dev.off()

###################################### Heatmap with Group labels - Transcripts ----
#Define the Experimental group colours
# Colours
# SW837_C9_CDS2_gRNA      SW837_C9_safe_gRNA     SW837_C9_uninfected     CDS2_Tumour_Cas9_WT CDS2_Tumour_Non-treated     CDS2_Tumour_Treated 
# "#A6CEE3"               "#559AC6"               "#3C8CAB"               "#94CA92"               "#7FC564"               "#33A02C" 
# "#f1b6da" "#c2a5cf"
# 
Exp_Groups=colorRampPalette(c(brewer.pal(n=4, name="Paired")))(length(unique(mdata$Group)))
names(Exp_Groups)<-unique(mdata$Group)
Cell_line=colorRampPalette(c("#f1b6da", "#c2a5cf"))(length(unique(mdata$Cell_line)))
names(Cell_line)<-unique(mdata$Cell_line)
group_cols<- list(Exp_Groups=Exp_Groups, 
                  Cell_line=Cell_line) # Set the group colors for the rest of the paper
#Defina the Heatmap annotation 
ha<- ComplexHeatmap::HeatmapAnnotation(Exp_Groups=mdata$Group[match(as.character(colnames(xc)), as.character(mdata$RNA_Sample_name))],
                                       Cell_line=mdata$Cell_line[match(as.character(colnames(xc)), as.character(mdata$RNA_Sample_name))],
                                       col = group_cols)
pdf(file.path(ov_fig_dir, paste("TPM_kallisto_transcripts_Pcor_heatmap_", studyid,"_samples_pGroups.pdf", sep="")),
    height = 12, width = 16)
ComplexHeatmap::Heatmap(t(xc), name = "Pearson_corr",
                        top_annotation = ha, 
                        col = hmcolors,
                        row_dend_width = unit(30, "mm"),
                        column_dend_height = unit(30, "mm"),
                        column_names_gp = grid::gpar(fontsize=10),
                        row_names_gp = grid::gpar(fontsize=10),
                        width =15, height = 10
)
dev.off()


#  Human Gene HEATMAP for the Pearson correlation comparison of all sample based on TPMs  ------------------
x<-tpm_tab
rownames(x)<- tpm_tab$target_id
x<-x[,2:dim(x)[2]]
#Calculate the LOG2
x<-log2(x+1)
#Calculate the pearson Correlation coefficient
xc<-cor(x)
hmcolors<-colorRampPalette(rev(brewer.pal(n=11,name = "RdYlBu" )))(70)

#GET the heatmap
pdf(file.path(ov_fig_dir, paste("TPM_kallisto_gene_Pcor_heatmap_", studyid,"_samples.pdf", sep="")),
    height = 12, width = 12)
gplots::heatmap.2(xc, col=hmcolors, main="Pearson correlation based on all\n genes \n log2(TPM+1)", 
                  key=TRUE, keysize=1.5,  cexCol=0.5, cexRow=0.5, trace="none", margins=c(7,7))
dev.off()
###################################### Heatmap with Group labels - Gene ----
#Define the Experimental group colours
# Colours
# SW837_C9_CDS2_gRNA      SW837_C9_safe_gRNA     SW837_C9_uninfected     CDS2_Tumour_Cas9_WT CDS2_Tumour_Non-treated     CDS2_Tumour_Treated 
# "#A6CEE3"               "#559AC6"               "#3C8CAB"               "#94CA92"               "#7FC564"               "#33A02C" 
Exp_Groups=colorRampPalette(c(brewer.pal(n=4, name="Paired")))(length(unique(mdata$Group)))
names(Exp_Groups)<-unique(mdata$Group)
Cell_line=colorRampPalette(c("#f1b6da", "#c2a5cf"))(length(unique(mdata$Cell_line)))
names(Cell_line)<-unique(mdata$Cell_line)
group_cols<- list(Exp_Groups=Exp_Groups, 
                  Cell_line=Cell_line) # Set the group colors for the rest of the paper
#Defina the Heatmap annotation 
ha<- ComplexHeatmap::HeatmapAnnotation(Exp_Groups=mdata$Group[match(as.character(colnames(xc)), as.character(mdata$RNA_Sample_name))], 
                                       Cell_line=mdata$Cell_line[match(as.character(colnames(xc)), as.character(mdata$RNA_Sample_name))],
                                       col = group_cols)
pdf(file.path(ov_fig_dir, paste("TPM_kallisto_gene_Pcor_heatmap_", studyid,"_samples_pGroups.pdf", sep="")),
    height = 12, width = 16)
ComplexHeatmap::Heatmap(t(xc), name = "Pearson_corr",
                        top_annotation = ha, 
                        col = hmcolors,
                        row_dend_width = unit(30, "mm"),
                        column_dend_height = unit(30, "mm"),
                        column_names_gp = grid::gpar(fontsize=10),
                        row_names_gp = grid::gpar(fontsize=10),
                        width =15, height = 10
)
dev.off()


# Function to calculate the GSEA for the Hallmark, C2 : CP and C6 cancer signatures 
get_hallmark_gsea<-function(so, de_gene_list, de_res_dir, prefix, sig_val=0.01, val_var="scaled_reads_per_base", gscategory="H"){
  if(!is.vector(de_gene_list)){ 
    stop("de_gene_list should be a vector")
  }
  if (gscategory=="H"){ # Select the 
    #Define GSEA database 
    pathwaysDF<- msigdbr::msigdbr("Homo sapiens", category = "H")
    pathways<-split(as.character(pathwaysDF$ensembl_gene), pathwaysDF$gs_name)
    gsea_file<-"Hallmark"
  } else if(gscategory=="C2"){
    #Define GSEA database 
    pathwaysDF<- msigdbr::msigdbr("Homo sapiens", category = "C2", subcategory = "CP")
    pathways<-split(as.character(pathwaysDF$ensembl_gene), pathwaysDF$gs_name)
    gsea_file<-"C2_CP"
  }else if(gscategory=="C6"){
    #Define GSEA database 
    pathwaysDF<- msigdbr::msigdbr("Homo sapiens", category = "C6")
    pathways<-split(as.character(pathwaysDF$ensembl_gene), pathwaysDF$gs_name)
    gsea_file<-"C6_Cancer_sig"
  } else{
    stop("Not a valid gscategory provided only H , C2 or C6 allowed")
  }
  
  # Val_var can be either "scaled_reads_per_base" or "tpm"
  # Do the GSEA analysis for the DE genes
  # Get the scaled values for GSEA
  tpm_tab_temp<-reshape2::dcast(so$obs_norm, target_id ~ sample, value.var = as.character(val_var))
  rownames(tpm_tab_temp)<-tpm_tab_temp$target_id
  tpm_tab_temp<-(tpm_tab_temp[, 2:dim(tpm_tab_temp)[2]])
  tpm_tab_temp<- tpm_tab_temp[rownames(tpm_tab_temp) %in% de_gene_list, ]
  tpm_tab_temp<-as.matrix(tpm_tab_temp)
  set.seed(as.numeric(Sys.getenv(envvar_seed, unset=default_seed)))
  gesecaRes <-geseca(pathways, tpm_tab_temp, minSize = 15, maxSize = 500, center=FALSE  )
  gesecaRes<-dplyr::filter(gesecaRes, padj <= sig_val)
  pheight<-8+round((dim(gesecaRes)[1]*10)/50,digits = 1)
  if(dim(gesecaRes)[1]==0){
    print("No enriched pathways found printing an empty plot", stdout())
    pdf(file = file.path(de_res_dir, paste(prefix, "_", gsea_file, "_GSEA_padj_", sig_val, "box.pdf", sep="")),
        height =pheight, width = 10 )
      ggplot() +
      geom_blank()
    dev.off()
  } else {
    gsplot<-plotGesecaTable(gesecaRes, pathways, E=tpm_tab_temp)
    ggsave(gsplot, file = file.path(de_res_dir, paste(prefix, "_", gsea_file,"_GSEA_padj_", sig_val, "box.pdf", sep="")),
           height =pheight, width = 10 )
  }
  data.table::fwrite(gesecaRes, 
                     file = file.path(de_res_dir, paste(prefix, "_", gsea_file, "_GSEA_padj_", sig_val, ".tsv", sep="")),
                     quote=FALSE, row.names = FALSE, sep="\t" )
}

######### DE with Sleuth for SW835 Cells
# perform_DE_pair function --- Development 
perform_DE_pair<-function(data_to_analyse, prefix, control, treatment, sig_val=0.01, de_res_dir, trans2gene, ncores ){
#  data_to_analyse  table for sleuth_prep
#  prefix<-"SW837_uninfected_VS_SW837_CDS2"  # Comparison name 
#  control<-"SW837_C9_uninfected" # Control name 
#  treatment<-"SW837_C9_CDS2_gRNA" # Group treatmente name
#  sig_val<-0.01  # Significance value 
#  de_res_dir<-file.path(results_dir, "DE_results", prefix) # Ouput dir for Comparison 
  if(!dir.exists(de_res_dir)){
    dir.create(de_res_dir, recursive = T)  
  }
  # To ensure Uninfected appeares first
  suppressWarnings(so <- sleuth_prep(data_to_analyse, ~CDS_groups, target_mapping = trans2gene, 
                    aggregation_column = "ens_gene", 
                    gene_mode=TRUE,
                    extra_bootstrap_summary=TRUE,
                    read_bootstrap_tpm = TRUE,
                    num_cores = ncores))
  # Fit the model with all the variables - using default
  so <- sleuth_fit(so, ~CDS_groups, "full" )
  # Fit the null model  - using default 
  so <- sleuth_fit(so, ~1, "reduced")
  # Likelihood ratio test
  so <- sleuth_lrt(so, "reduced", "full")
  so <- sleuth_wt(so, which_beta=paste0("CDS_groups", treatment)) 
  # Get the result table 
  sleuth_res_table<- sleuth_results(obj = so, test = "reduced:full", test_type = 'lrt', show_all = FALSE)
  sleuth_res_table <- dplyr::filter(sleuth_res_table, qval <= sig_val)
  
  # Save RDS object
  print("Generating Saving the sleuth object and results table plot ...", stdout())
  saveRDS(object = so, file =file.path(de_res_dir, paste0(prefix, "_sleuth_object.RDS")), compress = TRUE )
  #Save the results table 
  data.table::fwrite(sleuth_res_table, file = file.path(de_res_dir,paste0(prefix, "_DE_LRT_sig_", sig_val, ".tsv")), quote = FALSE, 
                                                       row.names=FALSE, col.names=TRUE, sep="\t")
  #PCA
  print("Generating the PCA plot ...", stdout())
  so_pca<-plot_pca(so, pc_x=1, pc_y=2, units = "scaled_reads_per_base", color_by = 'CDS_groups') +
    ggtitle(paste0(prefix,  " PCA"))
  so_pca_var<-plot_pc_variance(obj = so, pca_number = 2, units = "scaled_reads_per_base")+
    ggtitle(paste0(prefix,  " PCA var explained"))
  # Save PCAs
  pdf(file.path(de_res_dir, paste(prefix, "_PCA_scaledreads_per_base_varexp.pdf", sep="")),
      height = 12, width = 16)
  gridExtra::grid.arrange(so_pca, so_pca_var, ncol=2, nrow=1, widths=c(6,3))
  dev.off()
  
  # PLoT gene Jensen-Shannon divergence per samples
  print("Generating Jensen_Shannon_divergence plot ...", stdout())
  pdf(file.path(de_res_dir, paste(prefix, "_Jensen_Shannon_diverg_heatmap_scaledreads_per_base.pdf", sep="")),
      height = 8, width = 12)
  plot_sample_heatmap(obj = so)
  dev.off()
  
  return(so)
  
}

# Function to get the result table per Gene = models can be same as slueth allows and test_type should match
get_results_per_gene<-function(so, fname, sig_val=0.01, testv="reduced:full", test_typev = 'lrt', ens_hgnc=NULL ){
  if(is.null(ens_hgnc)){
    stop("Provide a table with the ensembl gene ids and hgnc")
  }
  # Get the result table 
  sleuth_resg_table<- sleuth_results(obj = so, test =testv , test_type = test_typev, show_all = FALSE)
  sleuth_resg_table <- dplyr::filter(sleuth_resg_table, qval <= sig_val)
  sleuth_res_table<-sleuth_resg_table
  sleuth_resg_table<-sleuth_resg_table[!duplicated(sleuth_resg_table$target_id),] 
  # Get the number of transcripts by counting number of times the gene ID appears on the table
  gene_list<-table(sleuth_res_table$target_id) 
  #Assign the number to the gene row
  sleuth_resg_table$num_diff_transcripts<- gene_list[sleuth_resg_table$target_id]
  #Add HGNC IDs
  sleuth_resg_table$hgnc_symbol<-ens_hgnc[match(sleuth_resg_table$target_id, ens_hgnc$ensembl_gene_id), c("hgnc_symbol")]
  sleuth_resg_table$hgnc_id<-ens_hgnc[match(sleuth_resg_table$target_id, ens_hgnc$ensembl_gene_id), c("hgnc_id")]
  sleuth_resg_table$chromosome_name<-ens_hgnc[match(sleuth_resg_table$target_id, ens_hgnc$ensembl_gene_id), c("chromosome_name")]
  sleuth_resg_table$start_position<-ens_hgnc[match(sleuth_resg_table$target_id, ens_hgnc$ensembl_gene_id), c("start_position")]
  sleuth_resg_table$end_position<-ens_hgnc[match(sleuth_resg_table$target_id, ens_hgnc$ensembl_gene_id), c("end_position")]
  #If test_type
  if(test_typev=="wt"){ # Log2FoldChange Calculation is Log2FC = b/ln(2), FoldChange=e^b
    #Use the beta values to get the log2Foldchange
    sleuth_resg_table$log2fC<- (sleuth_resg_table$b/log(2))
  }
  #Save the file where it got asked to 
  data.table::fwrite(sleuth_resg_table, file = fname, quote = FALSE, 
                     row.names=FALSE, col.names=TRUE, sep="\t")
  if(test_typev=="wt"){ # Print a log2 fold change filtered table
     temp_fname<-gsub(".tsv", "_log2FCFilt.tsv", fname, fixed = T)
     sleuth_resg_table<- sleuth_resg_table[(sleuth_resg_table$log2fC<=-0.5849625 | sleuth_resg_table$log2fC>=0.5849625 ),]
     data.table::fwrite(sleuth_resg_table[(sleuth_resg_table$log2fC<=-0.5849625 | sleuth_resg_table$log2fC>=0.5849625 ),], file = temp_fname, quote = FALSE, 
                        row.names=FALSE, col.names=TRUE, sep="\t")
  }
  return(sleuth_resg_table)
}


# DE SW837 Uninf 
comp_prefix<-"SW837_uninfected_VS_SW837_CDS2"
control_name<-"SW837_C9_uninfected"
treatment_name<-"SW837_C9_CDS2_gRNA"
fdrval<-0.01
de_res_dir_uninf_cds2<-file.path(results_dir, "DE_results", comp_prefix)
dir.create(de_res_dir_uninf_cds2, recursive = T)
data_to_analyse_input<-s2c_sw[ s2c_sw$CDS_groups %in% c( control_name, treatment_name), ]
data_to_analyse_input$CDS_groups<- factor(data_to_analyse_input$CDS_groups,levels = c(control_name, treatment_name))
nthreads=7

so_sw_uninf_cds2<- perform_DE_pair(data_to_analyse=data_to_analyse_input , 
                                   prefix=comp_prefix, 
                                   control=control_name,
                                   treatment=treatment_name,
                                   sig_val=fdrval,
                                   de_res_dir=file.path(results_dir, "DE_results", comp_prefix),
                                   trans2gene=t2g,
                                   ncores=nthreads )

temp_so<-so_sw_uninf_cds2
#Save the results table 
t1<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_LRT_sig_", fdrval, "_perGeneRes.tsv")),
                     sig_val = fdrval, testv="reduced:full", test_typev = 'lrt', ens_hgnc = ensv_genes)
dir.create(file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), recursive = TRUE)
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
t2<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_WaldTest_sig_", fdrval, "_perGeneRes.tsv")),
                     sig_val = fdrval, testv=paste0("CDS_groups", treatment_name), test_typev = 'wt',  ens_hgnc = ensv_genes)
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
# Plot  Distributions
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_scaledreads_per_base_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "scaled_reads_per_base",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_TPM_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "tpm",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
# Plot CDS2 scaled reads per base
print("Generating CDS2 plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "scaled_reads_per_base")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_FRRS1L_ENSG00000260230_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000260230", units = "tpm") + ggtitle("ENSG00000260230 FRRS1L")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SMG1_ENSG00000157106_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000157106", units = "tpm") + ggtitle("ENSG00000157106 SMG1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_CHAC1_DDIT4_SESN2_TPM_pergrp.pdf", sep="")),
    height = 20, width = 12)
a<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
b<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000128965", units = "tpm") + ggtitle("ENSG00000128965 CHAC1")
c<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000168209", units = "tpm")+ ggtitle("ENSG00000168209 DDIT4")
d<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000285069", units = "tpm") + ggtitle("ENSG00000285069 SESN2")
  gridExtra::grid.arrange(a, b,c,d,ncol=1, nrow=4)          
dev.off()
# QQ plots
print("Generating the QQ plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_QQ_plotPCA_scaledreads_per_base_lrt.pdf", sep="")),
    height = 12, width = 16)
plot_qq(temp_so, test = 'reduced:full', test_type = 'lrt', sig_level = fdrval)
dev.off()
# Plot MA with WT
print("Generating MA plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_MA_plot_WaldTest_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_ma(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
# Plot Volcano
print("Generating Volcano plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
png(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.png", sep="")),
    height = 400, width = 600)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval )
dev.off()
rm(temp_so, t1, t2, a,b,c,d)


# DE SW837 SW837_C9_safe_gRNA  
comp_prefix<-"SW837_C9_safe_gRNA_VS_SW837_CDS2"
control_name<-"SW837_C9_safe_gRNA"
treatment_name<-"SW837_C9_CDS2_gRNA"
fdrval<-0.01
de_res_dir_uninf_cds2<-file.path(results_dir, "DE_results", comp_prefix)
dir.create(de_res_dir_uninf_cds2, recursive = T)
data_to_analyse_input<-s2c_sw[ s2c_sw$CDS_groups %in% c( control_name, treatment_name), ]
data_to_analyse_input$CDS_groups<- factor(data_to_analyse_input$CDS_groups,levels = c(control_name, treatment_name))
nthreads=7

so_sw_safe_cds2<- perform_DE_pair(data_to_analyse=data_to_analyse_input , 
                                   prefix=comp_prefix, 
                                   control=control_name,
                                   treatment=treatment_name,
                                   sig_val=fdrval,
                                   de_res_dir=file.path(results_dir, "DE_results", comp_prefix),
                                   trans2gene=t2g,
                                   ncores=nthreads )
temp_so<-so_sw_safe_cds2
#Save the results table 
t1<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_LRT_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv="reduced:full", test_typev = 'lrt', ens_hgnc = ensv_genes)
dir.create(file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), recursive = TRUE)
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
t2<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_WaldTest_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv=paste0("CDS_groups", treatment_name), test_typev = 'wt',  ens_hgnc = ensv_genes)
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
# Plot  Distributions
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_scaledreads_per_base_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "scaled_reads_per_base",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_TPM_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "tpm",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
# Plot CDS2 scaled reads per base
print("Generating CDS2 plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "scaled_reads_per_base")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_PLD4_ENSG00000166428_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000166428", units = "tpm") + ggtitle("ENSG00000166428 PLD4")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_FRRS1L_ENSG00000260230_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000260230", units = "tpm") + ggtitle("ENSG00000260230 FRRS1L")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SMG1_ENSG00000157106_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000157106", units = "tpm") + ggtitle("ENSG00000157106 SMG1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_CHAC1_DDIT4_SESN2_TPM_pergrp.pdf", sep="")),
    height = 20, width = 12)
a<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
b<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000128965", units = "tpm") + ggtitle("ENSG00000128965 CHAC1")
c<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000168209", units = "tpm")+ ggtitle("ENSG00000168209 DDIT4")
d<-plot_bootstrap(obj = temp_so, target_id = "ENSG00000285069", units = "tpm") + ggtitle("ENSG00000285069 SESN2")
gridExtra::grid.arrange(a, b,c,d,ncol=1, nrow=4)          
dev.off()
# QQ plots
print("Generating the QQ plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_QQ_plotPCA_scaledreads_per_base_lrt.pdf", sep="")),
    height = 12, width = 16)
plot_qq(temp_so, test = 'reduced:full', test_type = 'lrt', sig_level = fdrval)
dev.off()
# Plot MA with WT
print("Generating MA plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_MA_plot_WaldTest_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_ma(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
# Plot Volcano
print("Generating Volcano plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
png(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.png", sep="")),
    height = 400, width = 600)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval )
dev.off()
rm(temp_so,t2, a,b,c,d)


# DE CDS2 Non_treated VS Treated 
comp_prefix<-"OMMCDS2_Tum_UnTreated_VS_CDS2_Tum_Treated"
control_name<-"CDS2_Tumour_Non-treated"
treatment_name<-"CDS2_Tumour_Treated"
failed_samples<-c("CDS2_Tumour_Cas9_WT_R4_RNA","CDS2_Tumour_R2_RNA", "CDS2_Tumour_Treated_R1_RNA")
fdrval<-0.01
de_res_dir_uninf_cds2<-file.path(results_dir, "DE_results", comp_prefix)
dir.create(de_res_dir_uninf_cds2, recursive = T)
data_to_analyse_input<-s2c_cds2[ s2c_cds2$CDS_groups %in% c( control_name, treatment_name), ]
# Removed failed samples for the expression analyses
data_to_analyse_input<- data_to_analyse_input[!(data_to_analyse_input$sample %in% failed_samples), ]
data_to_analyse_input$CDS_groups<- factor(data_to_analyse_input$CDS_groups,levels = c(control_name, treatment_name))
nthreads=7

so_om_untreat_cds2_treat<- perform_DE_pair(data_to_analyse=data_to_analyse_input , 
                                                     prefix=comp_prefix, 
                                                     control=control_name,
                                                     treatment=treatment_name,
                                                     sig_val=fdrval,
                                                     de_res_dir=file.path(results_dir, "DE_results", comp_prefix),
                                                     trans2gene=t2g,
                                                     ncores=nthreads )
temp_so<-so_om_untreat_cds2_treat
#Save the results table 
t1<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_LRT_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv="reduced:full", test_typev = 'lrt', ens_hgnc = ensv_genes)
# No GSEA as it has 0 DE genes found 
t2<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_WaldTest_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv=paste0("CDS_groups", treatment_name), test_typev = 'wt',  ens_hgnc = ensv_genes)
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
# Plot  Distributions
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_scaledreads_per_base_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "scaled_reads_per_base",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_TPM_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "tpm",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
# Plot CDS2 scaled reads per base
print("Generating CDS2 plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "scaled_reads_per_base")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_FRRS1L_ENSG00000260230_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000260230", units = "tpm") + ggtitle("ENSG00000260230 FRRS1L")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SMG1_ENSG00000157106_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000157106", units = "tpm") + ggtitle("ENSG00000157106 SMG1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS1_ENSG00000163624_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000163624", units = "tpm") + ggtitle("ENSG00000105392 CDS1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_OSGIN1_ENSG00000140961_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000140961", units = "tpm") + ggtitle("ENSG00000140961 OSGIN1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SNX19_ENSG00000120451_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000120451", units = "tpm") + ggtitle("ENSG00000120451 SNX19")
dev.off()

# QQ plots
print("Generating the QQ plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_QQ_plotPCA_scaledreads_per_base_lrt.pdf", sep="")),
    height = 12, width = 16)
plot_qq(temp_so, test = 'reduced:full', test_type = 'lrt', sig_level = fdrval)
dev.off()
# Plot MA with WT
print("Generating MA plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_MA_plot_WaldTest_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_ma(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
# Plot Volcano
print("Generating Volcano plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
png(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.png", sep="")),
    height = 400, width = 600)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval )
dev.off()
rm(temp_so, t2, t1)


# DE CDS2 Cas9 WT VS Treated 
comp_prefix<-"OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_Treated"
control_name<-"CDS2_Tumour_Cas9_WT"
treatment_name<-"CDS2_Tumour_Treated"
failed_samples<-c("CDS2_Tumour_Cas9_WT_R4_RNA","CDS2_Tumour_R2_RNA", "CDS2_Tumour_Treated_R1_RNA")
fdrval<-0.01
de_res_dir_uninf_cds2<-file.path(results_dir, "DE_results", comp_prefix)
dir.create(de_res_dir_uninf_cds2, recursive = T)
data_to_analyse_input<-s2c_cds2[ s2c_cds2$CDS_groups %in% c( control_name, treatment_name), ]
# Removed failed samples for the expression analyses
data_to_analyse_input<- data_to_analyse_input[!(data_to_analyse_input$sample %in% failed_samples), ]
data_to_analyse_input$CDS_groups<- factor(data_to_analyse_input$CDS_groups,levels = c(control_name, treatment_name))
nthreads=7

so_om_wt_cds2_treat<- perform_DE_pair(data_to_analyse=data_to_analyse_input , 
                                   prefix=comp_prefix, 
                                   control=control_name,
                                   treatment=treatment_name,
                                   sig_val=fdrval,
                                   de_res_dir=file.path(results_dir, "DE_results", comp_prefix),
                                   trans2gene=t2g,
                                   ncores=nthreads )
temp_so<-so_om_wt_cds2_treat
#Save the results table 
t1<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_LRT_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv="reduced:full", test_typev = 'lrt', ens_hgnc = ensv_genes)
dir.create(file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), recursive = TRUE)
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
t2<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_WaldTest_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv=paste0("CDS_groups", treatment_name), test_typev = 'wt',  ens_hgnc = ensv_genes)
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
# Plot  Distributions
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_scaledreads_per_base_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "scaled_reads_per_base",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_TPM_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "tpm",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
# Plot CDS2 scaled reads per base
print("Generating CDS2 plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "scaled_reads_per_base")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CRX_ENSG00000105392_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000105392", units = "tpm") + ggtitle("ENSG00000105392 CRX")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS1_ENSG00000163624_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000163624", units = "tpm") + ggtitle("ENSG00000105392 CDS1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_OSGIN1_ENSG00000140961_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000140961", units = "tpm") + ggtitle("ENSG00000140961 OSGIN1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SNX19_ENSG00000120451_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000120451", units = "tpm") + ggtitle("ENSG00000120451 SNX19")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_NRL_ENSG00000129535_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000129535", units = "tpm") + ggtitle("ENSG00000129535 NRL")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_FRRS1L_ENSG00000260230_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000260230", units = "tpm") + ggtitle("ENSG00000260230 FRRS1L")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SMG1_ENSG00000157106_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000157106", units = "tpm") + ggtitle("ENSG00000157106 SMG1")
dev.off()
# QQ plots
print("Generating the QQ plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_QQ_plotPCA_scaledreads_per_base_lrt.pdf", sep="")),
    height = 12, width = 16)
plot_qq(temp_so, test = 'reduced:full', test_type = 'lrt', sig_level = fdrval)
dev.off()
# Plot MA with WT
print("Generating MA plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_MA_plot_WaldTest_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_ma(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
# Plot Volcano
print("Generating Volcano plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
png(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.png", sep="")),
    height = 400, width = 600)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval )
dev.off()
rm(temp_so, t2, t1)


# DE CDS2 Cas9 WT VS UnTreated 
comp_prefix<-"OMMCDS2_Tum_Cas9WT_VS_CDS2_Tum_UnTreated"
control_name<-"CDS2_Tumour_Cas9_WT"
treatment_name<-"CDS2_Tumour_Non-treated"
failed_samples<-c("CDS2_Tumour_Cas9_WT_R4_RNA","CDS2_Tumour_R2_RNA", "CDS2_Tumour_Treated_R1_RNA")
fdrval<-0.01
de_res_dir_uninf_cds2<-file.path(results_dir, "DE_results", comp_prefix)
dir.create(de_res_dir_uninf_cds2, recursive = T)
data_to_analyse_input<-s2c_cds2[ s2c_cds2$CDS_groups %in% c( control_name, treatment_name), ]
# Removed failed samples for the expression analyses
data_to_analyse_input<- data_to_analyse_input[!(data_to_analyse_input$sample %in% failed_samples), ]
data_to_analyse_input$CDS_groups<- factor(data_to_analyse_input$CDS_groups,levels = c(control_name, treatment_name))
nthreads=7

so_om_wt_cds2_untreat<- perform_DE_pair(data_to_analyse=data_to_analyse_input , 
                                           prefix=comp_prefix, 
                                           control=control_name,
                                           treatment=treatment_name,
                                           sig_val=0.01,
                                           de_res_dir=file.path(results_dir, "DE_results", comp_prefix),
                                           trans2gene=t2g,
                                           ncores=nthreads )
temp_so<-so_om_wt_cds2_untreat
#Save the results table 
t1<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_LRT_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv="reduced:full", test_typev = 'lrt', ens_hgnc = ensv_genes)
dir.create(file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), recursive = TRUE)
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, de_res_dir=file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t1$target_id, file.path(de_res_dir_uninf_cds2, "Lrt_gsea"), 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
t2<-get_results_per_gene(so=temp_so, fname=file.path(de_res_dir_uninf_cds2,paste0(comp_prefix, "_DE_WaldTest_sig_", fdrval, "_perGeneRes.tsv")),
                         sig_val = fdrval, testv=paste0("CDS_groups", treatment_name), test_typev = 'wt',  ens_hgnc = ensv_genes)
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="H")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C2")
get_hallmark_gsea(so=temp_so, de_gene_list=t2$target_id, de_res_dir=de_res_dir_uninf_cds2, 
                  prefix=comp_prefix, sig_val=fdrval, val_var="scaled_reads_per_base", gscategory="C6")
# Plot  Distributions
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_scaledreads_per_base_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "scaled_reads_per_base",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_TPM_distribution.pdf", sep="")),
    height = 8, width = 9)
plot_group_density(temp_so, use_filtered = TRUE, units = "tpm",
                   trans = "log", grouping= setdiff(colnames(temp_so$sample_to_covariates),"sample"), offset = 1)
dev.off()
# Plot CDS2 scaled reads per base
print("Generating CDS2 plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_scaledreads_per_base_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "scaled_reads_per_base")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS2_ENSG00000101290_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000101290", units = "tpm") + ggtitle("ENSG00000101290 CDS2")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_FRRS1L_ENSG00000260230_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000260230", units = "tpm") + ggtitle("ENSG00000260230 FRRS1L")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SMG1_ENSG00000157106_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000157106", units = "tpm") + ggtitle("ENSG00000157106 SMG1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_CDS1_ENSG00000163624_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000163624", units = "tpm") + ggtitle("ENSG00000105392 CDS1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_OSGIN1_ENSG00000140961_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000140961", units = "tpm") + ggtitle("ENSG00000140961 OSGIN1")
dev.off()
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_SNX19_ENSG00000120451_TPM_pergrp.pdf", sep="")),
    height = 8, width = 12)
plot_bootstrap(obj = temp_so, target_id = "ENSG00000120451", units = "tpm") + ggtitle("ENSG00000120451 SNX19")
dev.off()

# QQ plots
print("Generating the QQ plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_QQ_plotPCA_scaledreads_per_base_lrt.pdf", sep="")),
    height = 12, width = 16)
plot_qq(temp_so, test = 'reduced:full', test_type = 'lrt', sig_level = fdrval)
dev.off()
# Plot MA with WT
print("Generating MA plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_MA_plot_WaldTest_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_ma(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
# Plot Volcano
print("Generating Volcano plot ...", stdout())
pdf(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.pdf", sep="")),
    height = 8, width = 12)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval)
dev.off()
png(file.path(de_res_dir_uninf_cds2, paste(comp_prefix, "_Volcano_treatment.png", sep="")),
    height = 400, width = 600)
plot_volcano(obj = temp_so, test=paste0("CDS_groups", treatment_name), sig_level=fdrval )
dev.off()
rm(temp_so, t2, t1)



# Save Rsession
save.image(file=file.path(DE_res_dir, "CDS2_targeting_kallisto_sleuth_env.RData" ), compress = TRUE )
