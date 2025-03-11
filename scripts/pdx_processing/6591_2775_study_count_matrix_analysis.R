############################################################
# 6709_2877_study_count_matrix_analysis_imp_from_canapps.R
# This script does the analysis for project 5823
#
# 0. Set the variables for the Project  
#   0.1 Merge the iRODS manifest with the metadata manifest to include the experimental information in a single file
#   1. DATA Loading
#     1.1 Create the master table with all the counts for the different samples of the project
#     1.2 Getting the ENSEMBL INFORMATION for retrieven ENSEMBLv104 data with BiomaRt
#       1.2.1 Getting the entire gene information form ENSEMBL INFORMATION for v104 
#       1.2.2 Getting the gene Length and  GC content  length information form ENSEMBL INFORMATION for all genes v104 
#       1.2.3 Adding the length of the ERCC spike ins to the table  
#   2. Define Functions for TPM calculations
#       get_TPM_from_master_table
#       get_nofeature
#       get_ambiguous
#       get_fcount_per_sample
#   3. Get the Transcript Per Million (TPMs) scaling
#     3.1 Genehumane the HEATMAP for the Pearson correlation comparison of all sample based on TPMs    
#     3.2 Save the table with the TPMs for all the samples with the respective gene_name and ENSEMBL gene ID (not all genes have gene names)   
#   4. Quality Checks 
#     4.1 Plot total number of counts per sample by category
#     4.2 Ploting total correlation of the total number of ambiguous reads and the RIN score
#     4.3 Remove samples that did  not pass QC (i.e. samples <40e6 reads counts)
# 5.0 Heatmap of Sample Pair pearson correlation based on from all expression in TPMS of all protein coding genes with the metadata of experimental information right next to it.
#     5.1 Fully annotated TPM heatmap with the metadata information right next to it. 
#     5.1 Heatmap with same information as above plus RIN QC information  
# 6.0 Compile Final table with STAR2.7.9a mapping results for all samples in the study
#   6.1 Mapping statistics plots ------------------
# Save the session  version Information  ------------------
#
#Created: 09/06/2022
#Author: Martin Del Castillo Velasco-Herrera - mdc1@sanger.ac.uk
###################################################################
# if you need to Install the  packages
if(!require("BiocManager", quietly = T)){ install.packages("BiocManager") }
bcpkgs<-c("ggplot2",
          "ggvenn",
          "RColorBrewer",
          "biomaRt",
          "tidyverse",
          "tidyr",
          "dplyr",
          "rmarkdown",
          "kableExtra",
          "DESeq2",
          "ComplexHeatmap",
          "gridExtra")
suppressMessages( BiocManager::install(setdiff(bcpkgs, rownames(installed.packages())), version = '3.14', update =FALSE) )
pkgs<-c( "data.table",
         "doParallel",
         "pheatmap",
         "shiny",
         "colourpicker",
         "DT",
         "plotly",
         "htmlwidgets",
         "patchwork",
         "png",
         "dendextend",
         "gplots",
         "cowplot")
suppressMessages( install.packages(setdiff(pkgs, rownames(installed.packages())), update =FALSE, repos ="https://cran.ma.imperial.ac.uk/") )

library(ggplot2)
library(RColorBrewer)
library(gplots)
library(biomaRt)
library(doParallel)
library(ggvenn)
library(DESeq2)
library(ComplexHeatmap)
library(reshape2)
library(cowplot)
library(tidyr)
library(dplyr)
library(data.table)
library(knitr)
library(R.utils)
library(patchwork)
library(limma)
#library(EDASeq)
#library(reshape2)



#############################################################
# 0. Set the variables for the Project  ------------------
print("setting the variables for the Project")
#FIRST the main project directory All other folders will be assumed to be located under this folder
#project_dir<-"/Users/mdc1/Desktop/team113sc124/projects/6591_PDX_models_Latin_America_RNAseq"
project_dir<-"/lustre/scratch124/casm/team113/projects/6591_PDX_models_Latin_America_RNAseq"

#Get the project name - Is the same as the project dir folder
project_name<-strsplit(project_dir, "/", )[[1]][length(strsplit(project_dir, "/", )[[1]])]
project_id<-strsplit(project_name, "_", )[[1]][1] # obtain the ID of the project
study_id<-"2775"

#Defining the directories for the analyses and inputs
# star_dir<-file.path(project_dir, "STAR2_bams")
htseqdir<-file.path(project_dir, "htseq_counts")
manif_dir<-file.path(project_dir, "manifests")
results_dir<-file.path(project_dir, "results")
res_qc_dir<- file.path(results_dir, "counts_qc")
xfiltres_dir<-file.path(results_dir, "xfilter_res")

#Create the results directory
dir.create(results_dir, recursive = T)
dir.create(res_qc_dir, recursive = T)
dir.create(xfiltres_dir, recursive = T)
setwd(results_dir)

#############################################################
# 0.1 Merge the iRODS manifest with the samples details from canapps ------------------
raw_manif<- as.data.frame(data.table::fread(file=file.path(manif_dir, "6591_cram_manifest_INFO_from_iRODS_wbam_counts_qc.txt"), sep="\t", header=TRUE))
raw_manif$seq_batch<- paste(raw_manif$id_run, raw_manif$lane, sep="_")
# dim(raw_manif)
# 
# #   0.1.1 Load the samples details from canapps sample details table ------------------
# #Load the sample details from CASM canapps for the project
# cp_sp_details<-read.csv(file = file.path(manif_dir, paste0(study_id, "_Cancer_Pipeline_Reports_SampleDetails.txt")), header = TRUE, sep="\t", stringsAsFactors = FALSE)
# 
# #   0.1.2 Load the samples details from canapps sample QC table ------------------
# #Load the sample details from CASM canapps  for the project
# cp_rnaqc_details<-read.csv(file = file.path(manif_dir, paste0(study_id,"_export-Automatically_generated_table_data.tsv")), header = TRUE, sep="\t", stringsAsFactors = FALSE)
# #Then we need to fix the names of the samples in the RNA QC because it mixes supplier sample name and the EGA STUDY sample name 
# cp_rnaqc_details$sample<-lapply(strsplit(cp_rnaqc_details$Sample, split =" ", fixed = T ), function(x){ x[[1]]} ) %>% unlist
# 
# #   0.1.3 Merge the samples details and QC table from canapps ------------------
# # Since the samples were sequenced using a single lane for each plex the merge doesn't need to take into account mulple lanes/runs for each sample
# #Then we merge the CGP tables removing the cgp_sp_details Seq.X column 
# ttab<-cp_sp_details[match(cp_rnaqc_details$sample, cp_sp_details$Sample), !(colnames(cp_sp_details) %in% c("Seq.X"))   ]
# #The merge removes the sample_name from the sample name  details(col1) from samp details table as this is already contained wihtin the QC table
# cp_master_tab<-cbind(cp_rnaqc_details , ttab[,c(2:dim(ttab)[2])] ) 
# rm(ttab)
# 
# #   0.1.4 Load the sequencescape information  ------------------
# sqp_manif<-read.csv(file = file.path(manif_dir, "5823_master_seqscape_manifest.txt"), header = TRUE, sep="\t", stringsAsFactors = FALSE)
# 
# #   0.1.4 Merge all the information CGP ans SEqscape into a single manifest ------------------
# #get the information form the sequencescape manifest to add into 
# ttab<-sqp_manif[ match(cp_master_tab$sample, sqp_manif$SUPPLIER_SAMPLE_NAME) ,]
# merged_cp_sqp_manif<- cbind(cp_master_tab, ttab[,c("SANGER_PLATE_ID","WELL","SANGER_SAMPLE_ID","SUPPLIER_SAMPLE_NAME","COHORT","GENDER", "seq_scape_manif_file","sample_sub_mision_batch")])
# 
# 
# # 0.2.Merge all the information into a manifest ------------------
# print("Merge the iRODS manifest with the merged metadata manifest to include the experimental information in a single file ")
merged_cp_sqp_manif<-read.csv(file = file.path(manif_dir, "6591_2775_master_manifest_wmetadata_for_proj_wNOD_bams_xfb_wcount.txt"), header = TRUE, sep="\t", stringsAsFactors = FALSE)
#This needs to be done onece the bam stats are collected
# master_manif<- cbind(raw_manif, merged_cp_sqp_manif[,])
master_manif<- merged_cp_sqp_manif
colnames(master_manif)[1]<-"Sanger_sample"

# write.table(master_manif,file=file.path(manif_dir,paste("Master_manifest_wmetadata_for_proj_",project_name, ".txt", sep="")), col.names = TRUE, row.names = FALSE, sep="\t")



#############################################################
# 1. DATA Loading ------------------
print("Count DATA loading")
#Reading the Sample Names
# master_manif<-raw_manif
# REad ONLY the names of samples that were sequencing completes 
sample_list<-master_manif$sample_supplier_name[grepl("Y", master_manif$Seq.Comp)]

#     1.1 Create the master table with all the counts for the different samples of the project ------------------
# Reading the tables  with the counts made with htseq using ENSEMBL v104
print("Create the counts master table ")
master_table<-NULL
#READ THE SAMPLES
print("Reading the count tables into a master matrix")
system.time(for(i in 1:length(sample_list)){
  temp<-NULL
  #Read the sample
  temp<-data.table::fread(file=file.path(htseqdir, paste(sample_list[i], "_htseq_counts_srverse.gz", sep="",collapse="")), header=FALSE, sep="\t", stringsAsFactors = F )
  temp<-as.data.frame(temp)
  #temp<-read.csv(file=file.path(htseqdir, paste(sample_list[i], "_htseq_counts_srverse.count.gz", sep="",collapse="")), header=FALSE, sep="\t", stringsAsFactors = F )
  #Then if this is the first sample we are reading we use the first column with the ENSEMBL gene IDs as the rowname and add the Sanger_ID as the column name
  if(i==1){
    master_table<-cbind(temp[,2])
    rownames(master_table)<-temp$V1 #Setting up the ENSEMBL Gene IDs as rownames and ERCCs
    colnames(master_table)<-sample_list[i]
  }else{
    master_table<-cbind(master_table, temp[,2])
    colnames(master_table)<-sample_list[1:i]
    rownames(master_table)<-temp$V1 #Setting up the ENSEMBL Gene IDs as rownames and ERCCs
  }
}
)

#Then we get the table for just the ERCC v92s spike ins in case they were used
ercc_mastertable<-master_table[grepl("ERCC", rownames(master_table)), ]

#Create the compiled count tables with the counts of the entire project
ENS_ID<-rownames(master_table) #temporary variable that will have the ENSEMBL IDS for the genes 
write.table(cbind(ENS_ID, master_table), file=file.path(results_dir,paste(project_name,"_HTSeq_Sreverse_STAR_ENSv103_ERCC_RAW_frag_counts_sample_supplier_name_IDs.tsv", sep="")), quote=FALSE, sep="\t", row.names = F)
#Create the count table with only the information for the ERCC spikes
ENS_ID<-rownames(ercc_mastertable) #temporary variable that will have the ENSEMBL IDS for the genes 
write.table(cbind(ENS_ID, ercc_mastertable), file=file.path(results_dir,paste(project_name,"_HTSeq_Sreverse_STAR_ERCC_ONLY_RAW_frag_counts_sample_sanger_IDs.tsv", sep="")), quote=FALSE, sep="\t", row.names = F)


#########################################################
#     1.2 Getting the gene Length information for the genes included for all genes used from ENSEMBL v103 ------------------
## Generally this is done with EDASeq package however, there is no way to fix the ENSEMBL version for EDASeq length
# so I'm using a file we generated for DERMATLAS from ENSEMBLv103
print("Getting the gene Length and  GC content  length information form ENSEMBL INFORMATION for all genes v103 ")
#Read current DERMATLAS gene metrics file used for TPM conversion
ens_gene_metrics<- read.csv(file=file.path(project_dir, "reference", "DERMATLAS-genes-GRCh38_103.txt"), header = T, sep="\t", stringsAsFactors = F)
colnames(ens_gene_metrics)[1]<- "ENSEMBL_GENE_ID"

#     1.3 Getting the ENSEMBL INFORMATION for the annotation used with BiomaRt ------------------
####################################### HERE
print("Import ENSEMBL v103 annotation information - if you are running this after ENSEMBv104 has been depracated change the human useMart line in the script")
#Get the name of all the gene IDs counted- this should be equal to the total num of rows minus (97) fields (which equate to the 92 ERCCs plus the last fields five
# added by htseq-counts which are  - "__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique"
ens_gene_ids<-rownames(master_table)[1:(length(row.names(master_table))-97) ]

#GET the genes.
# 07/09/2021  - ENSEMBL GENES 103 (SANGER UK)
#WHEN ENSEMBL changes to version 103 the way to lock the search for ENSEMBL v104 will be by uncommenting the command in the next line and
#running it instead of the command used originally
human<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="https://feb2021.archive.ensembl.org")
#     1.2.1 Getting the entire gene information form ENSEMBL INFORMATION for v103 ------------------
genebiotypes<-getBM(attributes=c("hgnc_symbol","hgnc_id","external_gene_name","external_gene_source", "ensembl_gene_id", "gene_biotype", "chromosome_name","start_position","end_position"), 
                    filters="ensembl_gene_id", values=ens_gene_ids,
                    mart=human, verbose = F)

#TEST TO CHECK that all the gene information downloaded is ok no repeated genes - NOTE : THER are two ens_gene_IDs that have two reported gene names for the same gene 
table(genebiotypes$ensembl_gene_id)[table(genebiotypes$ensembl_gene_id)>1]
# table(genebiotypes$ensembl_gene_id)[table(genebiotypes$ensembl_gene_id)>1]
# ENSG00000230417 ENSG00000254876 ENSG00000276085 
# 2               2               2 

#Get the names of the duplicates IDs (DUE to duplicated different gene names for the same ENSEMBL gene ID)
dupids<-names(table(genebiotypes$ensembl_gene_id)[table(genebiotypes$ensembl_gene_id)>1])
#Save a table of the  duplicated IDs
duptab<-genebiotypes[ genebiotypes$ensembl_gene_id %in% dupids,]
#To remove the duplicated entries and let the table have the same number of IDs in the count list
genebiotypes<-genebiotypes[!(c(1:length(genebiotypes$ensembl_gene_id)) %in% c(29396,42320, 53613 )), ]
table(rownames(master_table))[table(rownames(master_table))>1]

#Order the gene names in the the same order as the row names of the count tables 
genebiotypes<- genebiotypes[match(ens_gene_ids, genebiotypes$ensembl_gene_id),   ]
write.table(genebiotypes, file=file.path(results_dir,"ENSEMBL_v103_human_gene_biotypes_full_gene_information.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


#Get the list of protein coding genes
protein_cd_genes<-genebiotypes$ensembl_gene_id[genebiotypes$gene_biotype=="protein_coding"]
#Get the ribosomal genes only
rrna_cd_genes<-genebiotypes$ensembl_gene_id[genebiotypes$gene_biotype=="rRNA"]
#Logical vector to filterout all the protein coding genes only from master_table
cd_genes<-(rownames(master_table) %in% protein_cd_genes)
#Logical vector to get all rRNA genes
rrna_mtab_lf<-(rownames(master_table) %in% rrna_cd_genes)
#Get the information for all genes
all_genes<-genebiotypes


#For this I need to change EDASeq::getGeneLengthAndGCContent function to use the ensembl mart version i'm providing plus the new requirement of the specification of
#https:// section in the host URL

#To just call the default current ensembl use the function from the EDASeq::getGeneLengthAndGCContent 
#Since the Data contains a very large number of genes 60K and all the information and sequences are retrieved from ESNEMBL in a batch query. I split the
#data retrieval in 2 sections to make it manageable and prevent the function to break because it exceeds the allowed time with open connection for 
# The biomaRt access to ENSEMBL
# #set the 4 thresholds
# ens_gene_metrics<- NULL
# split_getGeneLengthAndGCContent<- function(ids,n=8, sp="hsa"){ # Function to split the number of IDs as input to getGeneLengthAndGCContent, 
#   # ids =ENSEMB_GENE IDS
#   # n= number of intervals to split the number of ids
#   # sp= species name in three letter used for biomaRt e.g "hsa" for Homo sapiens
#   out<-NULL
#   for (i in 1:n){
#   if(i==1){
#      tstart<-i
#      tend<-length(ens_gene_ids) %/% n
#      print(paste("Begining retrieval of section:", i, "of the IDs. Using IDs:",tstart," to ",tend, sep = " "))
#      out <- getGeneLengthAndGCContent(ids[tstart:tend], sp,mode="biomart" )
#      }else{
#        tstart<-tend+1
#        tend<-(length(ens_gene_ids)*i) %/% n
#        print(paste("Begining retrieval of section:", i, "of the IDs. Using IDs:", tstart, " to ",tend, sep=" "))
#        out<-rbind(out, getGeneLengthAndGCContent(ids[tstart:tend], sp,mode="biomart" ) )
#    } 
#   }
#   return(out)
# }
# 
# #Call my fucntion to break the lenght of the number of genes requested to ensemb by the EDASeq function 
# ens_gene_metrics<-split_getGeneLengthAndGCContent(ids=ens_gene_ids, n=8, sp="hsa")
# 
# #ens_gene_metrics<-getGeneLengthAndGCContent_m(ens_gene_ids, "hsa", mode = 'biomart')
# ens_gene_metrics<-as.data.frame(ens_gene_metrics)
# ens_gene_metrics$ENSEMBL_GENE_ID<- rownames(ens_gene_metrics)
# write.table(ens_gene_metrics, file=file.path(results_dir,"ENSEMBL_v103_human_genes_full_gene_length_GCcontent_information.txt"), quote=FALSE, row.names=F, col.names=TRUE, sep="\t")
# 
# #     1.2.3 Adding the length of the ERCC spike ins to the table  ------------------
# ercc_length<-read.csv(file = file.path(project_dir, "ERCC_Spike_in_info", "ERCCv92_gene_length.txt"), header = TRUE, sep="\t", stringsAsFactors = FALSE)
# row.names(ercc_length)<-ercc_length$ENSEMBL_GENE_ID
# ercc_length<-as.data.frame(ercc_length)
# ens_ercc_gene_length<- rbind(ens_gene_metrics[,c(3,1)],ercc_length[,1:2])
# write.table(ens_ercc_gene_length, file=file.path(results_dir,"ENSEMBL_v103_human_genes_ERCCv92_full_gene_length_information.txt"), quote=FALSE, row.names=F, col.names=TRUE, sep="\t")
# rm(ercc_length)


#     1.4 getting the Xenofilterstats table  ------------------
tpdx<- master_manif[master_manif$sample_type=="PDX", c("sample", "final_bam_used_counts")]
xfilt_stats<-NULL
i<-1
for(i in 1:dim(tpdx)[1]){
  tffile<-NULL
  temp<-NULL
  tffile<- gsub(".bam", ".bam.xfiltstats",tpdx[i,2], fixed = TRUE)
  temp<-read.csv(file = tffile, header = TRUE, sep = "\t", stringsAsFactors = F)
  xfilt_stats<- rbind(xfilt_stats, temp)
}
rm(temp, tffile, )
#Add metadata to the table
xfilt_stats<- cbind(xfilt_stats, master_manif[match(xfilt_stats$sample, master_manif$sample), c("sample_type", "tissue_phenotype")])
#Write table with the xenofilter stats for PDXs
write.table(xfilt_stats, file = file.path(xfiltres_dir, paste0(project_id, "_RNAseq_xenofilter_filtering_stats.tsv")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t" )


#######################################################################
# 2. Define Functions for TPM calculations   ------------------
#----------------------------------------------------------------------------
#   Create the function that will calculate the TPMs
#   TPM is very similar to RPKM and FPKM. The only difference is the order of opehumanions. Here’s how you calculate TPM:
#     1.-Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#     2.-Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#     3.-Divide the RPK values by the “per million” scaling factor. This gives you TPM.
#                   - --- IMPORTANT:  Row names should be ensembl Gene IDS ---
#--------------------------------------------------------------------------------
get_TPM_from_master_table<- function(tab=NULL, ense_gen_length=NULL, full_htseq_output=TRUE) {
  #Check if full_htseq_output is set (i.e. if the table contains the extra no-features and etc. )
  if(full_htseq_output){ 
    tpm<-tab[1:(dim(tab)[1]-5),]  
  }else {
    tpm<-tab
  }
  gne_len<-ense_gen_length
  #Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
  #Get the mean transcript lenght in KB for all the genes and in the same order
  genlengths<-gne_len$length[match(rownames(tpm), gne_len$ENSEMBL_GENE_ID)]
  genlengths<-(genlengths)/1000   
  #Divide the counts by transcript lenght in KB (RPK)
  for(i in 1:dim(tpm)[2]){
    tpm[,i]<-(tpm[,i]/genlengths) #get the reads per Kilobase
    #Get the number of RPK per samples
    tmp_rpk<- sum(tpm[,i])
    #Scale that per million
    tmp_rpk<- tmp_rpk/1000000
    #Divide the RPK values by the per million scaling factor
    tpm[,i]<-(tpm[,i]/tmp_rpk) #get the reads per Kilobase
  }
  return(tpm)
}
#Total number of no_feature or Ambgious:
#Function to get the No feature reads=
get_nofeature<-function(x){
  temp<-x[rownames(x) %in% c("__no_feature" ),]
  return(temp)
}
#Ambigious
get_ambiguous<-function(x){
  temp<-x[rownames(x) %in% c("__ambiguous" ),]
  return(temp)
}
#Alignment not_unique
get_alignment_not_unique<-function(x){
  temp<-x[rownames(x) %in% c("__alignment_not_unique" ),]
  return(temp)
}
#Total counts from htseq matrix
get_fcount_per_sample<-function(x){
  temp<-colSums(x[1:(dim(x)[1]-5),])
  return(temp)
}

get_rRNA_couns_persample<-function(x,rrna_gids){
  temp<-x[1:(dim(x)[1]-5),]
  temp<-colSums(temp[(rownames(temp) %in% rrna_gids),])
  return(temp)
}

##################################################################
# 3. Get the TPMs   ------------------
########################################################################################
print("Calculating the TPMs values for all the samples using the function withing this script ")
# #Calculate the TPMS
all_samples_TPMs<- get_TPM_from_master_table(tab = master_table, ense_gen_length =ens_gene_metrics, full_htseq_output = TRUE )

# Check that the TPMs from CASM sum all the columns
colSums(all_samples_TPMs)

# 
# #    3.1 Genehumane the HEATMAP for the Pearson correlation comparison of all sample based on TPMs  ------------------
x<-all_samples_TPMs
#Calculate the LOG2
x<-log2(x+1)
#Calculate the pearson Correlation coefficient
xc<-cor(x)
hmcolors<-colorRampPalette(rev(brewer.pal(n=11,name = "RdYlBu" )))(60)

#GET the heatmap
pdf(file.path(results_dir, paste("TPM_M_Pcor_heatmap_HTSeq_counts_STAR_", project_id,"_samples.pdf", sep="")))
heatmap.2(xc, col=hmcolors, main="Pearson correlation based on all\n genes\n log2(TPM+1)", 
          key=TRUE, keysize=1.25,  cexCol=0.55, cexRow=0.55, trace="none", margins=c(5,5))
dev.off()

#     3.2 Save the table with the TPMs for all the samples with the respective gene_name and ENSEMBL gene ID (not all genes have gene names)   ------------------
#Write the TPM table with two columns before the TPM values, gene_name, then ENSEMBL_gene_ID, then TPM values
tpmtab<-cbind( genebiotypes$external_gene_name[match(rownames(all_samples_TPMs), genebiotypes$ensembl_gene_id)], rownames(all_samples_TPMs), all_samples_TPMs)
colnames(tpmtab)[1:2]<-c("external_gene_name","ENSEMBL_GENE_ID")
write.table(tpmtab,
            file=file.path(results_dir,paste(project_id,"_HTSeq_STAR_ENSv103_ERCC_TPM_M_ENS_hgnc_symbol_IDS_ALLSAMPLES.txt")), quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
rm(tpmtab)


##########################################################
# 4. Quality Checks on samples    ------------------
print("Performing quality checks")
#Get the count types as a QC for  samples
ambiguous_reads<- get_ambiguous(master_table)
no_feature_reads<- get_nofeature(master_table)
gene_read_pair_counts<- get_fcount_per_sample(master_table)
rRNA_counts<-get_rRNA_couns_persample(master_table, rrna_gids = rrna_cd_genes)
alignment_not_unique<-get_alignment_not_unique(master_table)
rRNA_percent_counts<- (rRNA_counts/gene_read_pair_counts)*100

#Fail samples if the number of ambiguous reads is bigger than the number of gene pair read counts
fail_status<-NULL
qc_status_sum<-NULL
fail_status<-ifelse((gene_read_pair_counts<(ambiguous_reads+no_feature_reads)), "*", NA) # vector to show failed samples in plot
qc_status_sum<-ifelse((gene_read_pair_counts<(ambiguous_reads+no_feature_reads)), "high_ambig_readc", "pass") #Variable to keep a summary of reasons for failing QC

#Fail samples if they had less than 20Million reads counted
# For RNAseq expresison assessment sequencing depth is an important factor to be able to assess expression reliably. 
# We recommend to exclude samples with low number of read pairs (fragments) counted **<20 Million reads counted**. This would allow a reliable assessment of highly expressed genes (see **[Tarazona, S. _et. al._ 2011](https://genome.cshlp.org/content/21/12/2213.full) **).   

fail_status<- ifelse((fail_status %in% "*")|ifelse(gene_read_pair_counts<20e6,TRUE, FALSE) , #This condition is to check that no previously failed samples get overwritten  to NA
                     yes = "*" , no = NA )
for(i in 1:length(gene_read_pair_counts)){  # Code to update the qc_satus_sum 
  if(gene_read_pair_counts[i]<20e6){
    if(qc_status_sum[i]=="pass"){
      qc_status_sum[i]<-"low_rc"
    }else{
      qc_status_sum[i]<- paste(qc_status_sum[i], "low_rc", sep=",")
    }
  }
}


#plot the  Number of counted fragmens
count_table_type_sum<-cbind(names(gene_read_pair_counts), gene_read_pair_counts, no_feature_reads, ambiguous_reads,  alignment_not_unique, qc_status_sum,rRNA_counts,rRNA_percent_counts)
count_table_type_sum<-as.data.frame(count_table_type_sum)
colnames(count_table_type_sum)[1]<-"sample"
count_table_type_sum$gene_read_pair_counts<- as.numeric(as.character(count_table_type_sum$gene_read_pair_counts))
count_table_type_sum$no_feature_reads<- as.numeric(as.character(count_table_type_sum$no_feature_reads))
count_table_type_sum$ambiguous_reads<- as.numeric(as.character(count_table_type_sum$ambiguous_reads))
count_table_type_sum$alignment_not_unique<- as.numeric(as.character(count_table_type_sum$alignment_not_unique))

#Change to long format for ggplots
counts_sum_type<-reshape2::melt(data = count_table_type_sum[,1:(dim(count_table_type_sum)[2]-3)], id.vars ="sample")

#   4.1 Plot total number of counts per sample by category  ------------------
#Plot the total number of counts per sample  by category
count_sum_bplot <-ggplot(data=counts_sum_type, aes(x=sample, y=value ))+
  geom_bar(stat = "identity", aes(fill=variable), position = "dodge")+
  theme_bw()+
  scale_fill_brewer(type="qual",palette = 3  )+
  ggtitle(paste(project_name, "\nRead Pairs counted by HTseq by category"))+
  ylab("Counts")+
  xlab("Sample")+
  theme( axis.text.x = element_text(angle = 85 ,vjust =0.4, size = 10), 
         axis.title.x = element_text(face="bold", size = 12),
         axis.title.y = element_text(face="bold", size = 12),
         plot.title = element_text(hjust = 0.5, size = 14),
         axis.text.y = element_text( size = 12))+
  guides(fill=guide_legend(title="HTseq count types"))+
  geom_text(data=count_table_type_sum, aes(x=sample, y=gene_read_pair_counts,label=fail_status),angle=45, size=8,col="red")#To add the label of the Fail samples 
#count_sum_bplot
ggsave(count_sum_bplot, filename =file.path(res_qc_dir,paste(project_id, "_data_HTseq_count_summary_by_category_per_sample.pdf", sep="")),dpi=300, width = 20, height = 8)


#     4.3 Remove samples that did  not pass QC  ------------------
#Samples to Remove, FAILED samples (samples wiht higher no_feature or ambigious counts than uniquely mapped counts ) but in this case no samples had that 
#So this filter was removed
#Get the position of the samples that failed QC 
failedqc_samples<- !(count_table_type_sum$qc_status_sum %in% "pass")

#PLot the QC failing reasons 
t_qcsum<-cbind(names(table(count_table_type_sum$qc_status_sum)),table(count_table_type_sum$qc_status_sum)) 
colnames(t_qcsum)<- c("qc_status","value")
t_qcsum<- as.data.frame(t_qcsum,row.names =1:dim(t_qcsum)[1] ) 
t_qcsum$value<- as.numeric(t_qcsum$value)
qc_sump<- ggplot(data=t_qcsum, aes(x=qc_status, y=value ))+
  geom_bar(stat = "identity", aes(fill=qc_status), position = "dodge")+
  theme_bw()+
  scale_fill_brewer(type="qual",palette = "Set1"  )+
  ggtitle(paste(project_name, "\n QC summary"))+
  ylab("Number of Samplest")+
  theme( axis.text.x = element_blank(), 
         axis.title.x = element_text(face="bold", size = 14),
         axis.title.y = element_text(face="bold", size = 14),
         plot.title = element_text(hjust = 0.5),
         axis.text.y = element_text( size = 13))+
  guides(fill=guide_legend(title="QC status"))
# qc_sump
ggsave(qc_sump, filename =file.path(res_qc_dir,paste(project_id, "_bplot_QCsummary_by_qcstatus.pdf", sep="")),dpi=200, width = 7, height =8)
rm(t_qcsum, qc_sump)


#Save the table with the number of extra counts 
#Add the sample type information 
count_table_type_sum$sample_type<- master_manif$sample_type[match(count_table_type_sum$sample, master_manif$sample) ]
count_table_type_sum$Passed_QC<- FALSE
count_table_type_sum$Passed_QC[grepl("pass", count_table_type_sum$qc_status_sum, fixed = T)]<- TRUE
write.table(count_table_type_sum, file=file.path(res_qc_dir, paste(project_id ,"HTseq_counts_summary_by_class_and_QCsum.tsv", sep="")), col.names = TRUE, row.names=FALSE, sep="\t", quote = FALSE)

#get the final samples by removing the ones that FAILED QC of total number of reads in genes
master_table_final<- master_table[1:(dim(master_table)[1]-5),]
master_table_final<- master_table[, colnames(master_table)[(!failedqc_samples)]]
master_tpm_cp_final<- all_samples_TPMs[,colnames(all_samples_TPMs)[(!failedqc_samples)] ]  #For the TPMs calculated here


#     4.4 Save the table with the FINAL QC Passed M TPMs for all the samples with the respective gene_name and ENSEMBL gene ID (not all genes have gene names)   ------------------
#Write the TPM table with two columns before the TPM values, gene_name, then ENSEMBL_gene_ID, then TPM values
tpmtab<-cbind( genebiotypes$external_gene_name[match(rownames(master_tpm_cp_final), genebiotypes$ensembl_gene_id)], rownames(master_tpm_cp_final), master_tpm_cp_final)
colnames(tpmtab)[1:2]<-c("external_gene_name","ENSEMBL_GENE_ID")
write.table(tpmtab,
            file=file.path(results_dir,paste(project_id,"_M_TPMS_HTSeq_STAR_ENSv103_ERCC_TPM_M_ENS_hgnc_symbol_IDS_QC_PASSED.txt")), quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
rm(tpmtab)

#     4.6 Save the table with the FINAL QC Passed counts table  for all the samples with the respective gene_name and ENSEMBL gene ID (not all genes have gene names)   ------------------
#Write the count table with two columns before the count values, gene_name, then ENSEMBL_gene_ID, then TPM values
tpmtab<-cbind( genebiotypes$external_gene_name[match(rownames(master_table_final), genebiotypes$ensembl_gene_id)], rownames(master_table_final), master_table_final)
colnames(tpmtab)[1:2]<-c("external_gene_name","ENSEMBL_GENE_ID")
write.table(tpmtab,
            file=file.path(results_dir,paste(project_id,"_Counts_HTSeq_STAR_ENSv103_ERCC_TPM_M_ENS_hgnc_symbol_IDS_QC_PASSED.txt")), quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
rm(tpmtab)

#     4.7 ADD to the Master_manif the QC status and save   ------------------
Passed_QC<-vector(mode = "logical", length(master_manif$sample))
## Mark the samples that passed the QC  on the count_Table summary 
Passed_QC[ match(count_table_type_sum$sample[count_table_type_sum$qc_status_sum=="pass"], master_manif$sample) ]<- TRUE
master_manif$Passed_QC<- Passed_QC
write.table(master_manif, file = file.path(manif_dir, "6591_2775_master_manifest_wmetadata_for_proj_wNOD_bams_xfb_wcount_countQCPass.tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

##########################################################
# 5.0 Heatmap of Sample Pair pearson correlation based on from all TPMS across all gene with the metadata of experimental information right next to it. ------------------
print("Getting Correlation heatmaps")

#READ the metadata about tissue type and add it to the count table 
ttype<- read.csv(file= file.path(manif_dir,"cgp_contact", "6591_cgp_contact_master_manifest_ids_ony_samp_type.tsv"), header = T, stringsAsFactors = F, sep="\t")
ttype<- ttype[(ttype$human_sample_type=="RNA"), ]
count_table_type_sum$Histology<- ttype$tissue_histology[match(count_table_type_sum$sample, ttype$PD_ID)  ]

pdf(file.path(results_dir, paste("TPM_M_Pcor_heatmap_HTSeq_counts_STAR_", project_id,"_samples_Exp_group_annotation_QCpass.pdf", sep="")), width = 12, height = 10)
#Filter just protein coding genes
# x<-all_samples_TPMs[ match(protein_cd_genes, row.names(all_samples_TPMs)),]
x<-all_samples_TPMs
#Calculate the LOG2 (ALREAdy DONE that above)
x<-log2(x+1)
#Calculate the pearson Correlation coefficient
xc<-cor(x)

#Define the correlation colour spectrum
hmcolors<-colorRampPalette(rev(brewer.pal(n=11,name = "RdYlBu" )))(60)

#Define the sample type  group colours
Sample_type=c(brewer.pal(n=12, name="Paired")[c(8,10)])
names(Sample_type)<-unique(count_table_type_sum$sample_type)
Passed_QC=c(brewer.pal(n=12, name="Paired")[c(4,6)] )
names(Passed_QC)<-unique(count_table_type_sum$Passed_QC)
Histology=c(c(brewer.pal(n=12, name="Set3")[1:length(unique(count_table_type_sum$Histology))] ))
names(Histology)<-unique(count_table_type_sum$Histology)
group_cols<- list(Sample_type=Sample_type, 
                  Passed_QC=Passed_QC, 
                  Histology=Histology)

#Defina the Heatmap annotation 
tempv<-count_table_type_sum$sample_type
names(tempv)<- count_table_type_sum$sample
tempqc<-count_table_type_sum$Passed_QC
names(tempqc)<- count_table_type_sum$sample
temph<-count_table_type_sum$Histology
names(temph)<- count_table_type_sum$sample
ha<- HeatmapAnnotation(Sample_type=tempv,
                       Histology=temph, 
                       Passed_QC=tempqc,
                       col = group_cols)

Heatmap(t(xc),
        top_annotation = ha,
        name= "Pearson_corr",
        col = hmcolors,
        column_title = paste0(project_id, " PDX RNAseq\n","Pearson correlation based on all\n genes log2(TPM+1)"),
        row_dend_width = unit(30, "mm"),
        column_dend_height = unit(30, "mm"),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 9)),
        column_names_gp = gpar(fontsize=10),
        row_names_gp = gpar(fontsize=10),
        width =15, height = 10
)
dev.off()

##########################################################
# 5.1 Heatmap of Sample Pair pearson correlation based on from all TPMS across all gene with the metadata of experimental information right next to it. ------------------
print("Getting Correlation heatmaps")

pdf(file.path(results_dir, paste("TPM_M_Pcor_heatmap_HTSeq_counts_STAR_", project_id,"_samples_Prot_genes_Exp_group_annotation_QCpass.pdf", sep="")), width = 12, height = 10)
#Filter just protein coding genes
x<-all_samples_TPMs[ match(protein_cd_genes, row.names(all_samples_TPMs)),]
# x<-all_samples_TPMs
#Calculate the LOG2 (ALREAdy DONE that above)
x<-log2(x+1)
#Calculate the pearson Correlation coefficient
xc<-cor(x)

#Define the correlation colour spectrum
hmcolors<-colorRampPalette(rev(brewer.pal(n=11,name = "RdYlBu" )))(60)

#Define the sample type  group colours
Sample_type=c(brewer.pal(n=12, name="Paired")[c(8,10)])
names(Sample_type)<-unique(count_table_type_sum$sample_type)
Passed_QC=c(brewer.pal(n=12, name="Paired")[c(4,6)] )
names(Passed_QC)<-unique(count_table_type_sum$Passed_QC)
Histology=c(c(brewer.pal(n=12, name="Set3")[1:length(unique(count_table_type_sum$Histology))] ))
names(Histology)<-unique(count_table_type_sum$Histology)
group_cols<- list(Sample_type=Sample_type, 
                  Passed_QC=Passed_QC, 
                  Histology=Histology)

#Defina the Heatmap annotation 
tempv<-count_table_type_sum$sample_type
names(tempv)<- count_table_type_sum$sample
tempqc<-count_table_type_sum$Passed_QC
names(tempqc)<- count_table_type_sum$sample
temph<-count_table_type_sum$Histology
names(temph)<- count_table_type_sum$sample
ha<- HeatmapAnnotation(Sample_type=tempv,
                       Histology=temph, 
                       Passed_QC=tempqc,
                       col = group_cols)


Heatmap(t(xc),
        top_annotation = ha,
        name= "Pearson_corr",
        col = hmcolors,
        column_title = paste0(project_id, " PDX RNAseq\n","Pearson correlation based on all\n prot_coding genes log2(TPM+1)"),
        row_dend_width = unit(30, "mm"),
        column_dend_height = unit(30, "mm"),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 9)),
        column_names_gp = gpar(fontsize=10),
        row_names_gp = gpar(fontsize=10),
        width =15, height = 10
)
dev.off()

##########################################################
# 6.0 Heatmap of Sample Pair pearson correlation based on from all TPMS across all gene with the metadata of experimental information right next to it. ------------------
print("Getting Correlation heatmaps")

pdf(file.path(results_dir, paste("TPM_M_Pcor_heatmap_HTSeq_counts_STAR_", project_id,"_samples_Exp_group_annotation_QCpass_only.pdf", sep="")), width = 12, height = 10)
#Filter just protein coding genes
x<-all_samples_TPMs[, count_table_type_sum$sample[count_table_type_sum$Passed_QC] ]
#Calculate the LOG2 (ALREAdy DONE that above)
x<-log2(x+1)
#Calculate the pearson Correlation coefficient
xc<-cor(x)

#Define the correlation colour spectrum
hmcolors<-colorRampPalette(rev(brewer.pal(n=11,name = "RdYlBu" )))(60)

#Define the sample type  group colours
Sample_type=c(brewer.pal(n=12, name="Paired")[c(8,10)])
names(Sample_type)<-unique(count_table_type_sum$sample_type)
Passed_QC=c(brewer.pal(n=12, name="Paired")[c(4,6)] )
names(Passed_QC)<-unique(count_table_type_sum$Passed_QC)
Histology=c(c(brewer.pal(n=12, name="Set3")[1:length(unique(count_table_type_sum$Histology))] ))
names(Histology)<-unique(count_table_type_sum$Histology)
group_cols<- list(Sample_type=Sample_type, 
                  Passed_QC=Passed_QC, 
                  Histology=Histology)

#Defina the Heatmap annotation 
tempv<-count_table_type_sum$sample_type[count_table_type_sum$Passed_QC]
names(tempv)<- count_table_type_sum$sample[count_table_type_sum$Passed_QC]
tempqc<-count_table_type_sum$Passed_QC[count_table_type_sum$Passed_QC]
names(tempqc)<- count_table_type_sum$sample[count_table_type_sum$Passed_QC]
temph<-count_table_type_sum$Histology[count_table_type_sum$Passed_QC]
names(temph)<- count_table_type_sum$sample[count_table_type_sum$Passed_QC]
ha<- HeatmapAnnotation(Sample_type=tempv,
                       Histology=temph, 
                       Passed_QC=tempqc,
                       col = group_cols)
Heatmap(t(xc),
        top_annotation = ha,
        name= "Pearson_corr",
        col = hmcolors,
        column_title = paste0(project_id, " PDX RNAseq Pass_only \n","Pearson correlation based on all\n genes log2(TPM+1)"),
        row_dend_width = unit(30, "mm"),
        column_dend_height = unit(30, "mm"),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 9)),
        column_names_gp = gpar(fontsize=10),
        row_names_gp = gpar(fontsize=10),
        width =15, height = 10
)
dev.off()

##########################################################
# 6.1 Heatmap of Sample Pair pearson correlation based on from all TPMS across all gene with the metadata of experimental information right next to it. ------------------
print("Getting Correlation heatmaps")

pdf(file.path(results_dir, paste("TPM_M_Pcor_heatmap_HTSeq_counts_STAR_", project_id,"_samples_Prot_genes_Exp_group_annotation_QCpass_only.pdf", sep="")), width = 12, height = 10)
#Filter just protein coding genes
x<-all_samples_TPMs[ match(protein_cd_genes, row.names(all_samples_TPMs)),  count_table_type_sum$sample[count_table_type_sum$Passed_QC] ]
# x<-all_samples_TPMs
#Calculate the LOG2 (ALREAdy DONE that above)
x<-log2(x+1)
#Calculate the pearson Correlation coefficient
xc<-cor(x)

#Define the correlation colour spectrum
hmcolors<-colorRampPalette(rev(brewer.pal(n=11,name = "RdYlBu" )))(60)

#Define the sample type  group colours
Sample_type=c(brewer.pal(n=12, name="Paired")[c(8,10)])
names(Sample_type)<-unique(count_table_type_sum$sample_type)
Passed_QC=c(brewer.pal(n=12, name="Paired")[c(4,6)] )
names(Passed_QC)<-unique(count_table_type_sum$Passed_QC)
Histology=c(c(brewer.pal(n=12, name="Set3")[1:length(unique(count_table_type_sum$Histology))] ))
names(Histology)<-unique(count_table_type_sum$Histology)
group_cols<- list(Sample_type=Sample_type, 
                  Passed_QC=Passed_QC, 
                  Histology=Histology)

#Defina the Heatmap annotation 
tempv<-count_table_type_sum$sample_type[count_table_type_sum$Passed_QC]
names(tempv)<- count_table_type_sum$sample[count_table_type_sum$Passed_QC]
tempqc<-count_table_type_sum$Passed_QC[count_table_type_sum$Passed_QC]
names(tempqc)<- count_table_type_sum$sample[count_table_type_sum$Passed_QC]
temph<-count_table_type_sum$Histology[count_table_type_sum$Passed_QC]
names(temph)<- count_table_type_sum$sample[count_table_type_sum$Passed_QC]
ha<- HeatmapAnnotation(Sample_type=tempv,
                       Histology=temph, 
                       Passed_QC=tempqc,
                       col = group_cols)


Heatmap(t(xc),
        top_annotation = ha,
        name= "Pearson_corr",
        col = hmcolors,
        column_title = paste0(project_id, " PDX RNAseq\n","Pearson correlation based on all\n prot_coding genes log2(TPM+1)"),
        row_dend_width = unit(30, "mm"),
        column_dend_height = unit(30, "mm"),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 9)),
        column_names_gp = gpar(fontsize=10),
        row_names_gp = gpar(fontsize=10),
        width =15, height = 10
)
dev.off()




#######################################################################################################################
#Save the session  version INformation  ------------------
writeLines(capture.output(sessionInfo()), file.path(results_dir,"Analysis_R_session_INFO.txt"))



