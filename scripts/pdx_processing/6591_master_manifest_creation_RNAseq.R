############################################################
# 6591_master_manifest_creation_RNAseq.R
#
#Created: 01-07-2022
#Author: Martin Del Castillo Velasco-Herrera - mdc1@sanger.ac.uk
###################################################################
library(tidyr)
library(dplyr)

#############################################################
# 0. Set the variables for the Project  ------------------
#FIRST the main project directory All other folders will be assumed to be located under this folder
project_dir<-"/Users/mdc1/Documents/Projects/PDX_models_Latin_America/6591_PDX_models_from_Latin_America_RNAseq"
# project_dir<-"/Users/mdc1/Desktop/team113sc124/projects/6709_2877_DERMATLAS_Melanoma_MPNST_Cleveland_RNAseq"
#Get the project name - Is the same as the project dir folder
project_name<-strsplit(project_dir, "/", )[[1]][length(strsplit(project_dir, "/", )[[1]])]
project_id<-strsplit(project_name, "_", )[[1]][1] # obtain the ID of the project
study_id<-"2775"

#Defining the directories for the analyses and inputs
manif_dir<-file.path(project_dir, "manifests")
results_dir<-file.path(project_dir, "results")
cgp_contact_manifdir<-file.path(manif_dir, "cgp_contact")
seqscape_manifdir<-file.path(manif_dir, "sequencescape")

#Create the results directory
dir.create(results_dir, recursive = T)
setwd(results_dir)

#############################################################
# Load all of the required files
# 0.1 Merge the iRODS manifest with the samples details from canapps ------------------
raw_manif<- as.data.frame(data.table::fread(file=file.path(manif_dir, paste0(project_id,"_cram_manifest_INFO_from_iRODS.txt")), sep="\t", header=TRUE))
raw_manif$seq_batch<- paste(raw_manif$id_run, raw_manif$lane, sep="_")
dim(raw_manif)

#   0.1.1 Load the samples details from canapps sample details table ------------------
#Load the sample details from CASM canapps for the project
cp_sp_details<-read.csv(file = file.path(manif_dir, paste0(study_id, "_Cancer_Pipeline_Reports_SampleDetails.txt")), header = TRUE, sep="\t", stringsAsFactors = FALSE)
dim(cp_sp_details)
#   0.1.2 Load the samples details from canapps sample QC table ------------------
#Load the sample details from CASM canapps  for the project
cp_qc_details<-read.csv(file = file.path(manif_dir, paste0(study_id,"_export-Automatically_generated_table_data.csv")), header = TRUE, stringsAsFactors = FALSE)
#Then we need to fix the names of the samples in the RNA QC because it mixes supplier sample name and the EGA STUDY sample name 
cp_qc_details$sample<-lapply(strsplit(cp_qc_details$Sample, split =" ", fixed = T ), function(x){ x[[1]]} ) %>% unlist
dim(cp_qc_details)
#Load information about the samples whether they are PDX or not
metadata_tab<-read.table(file = file.path(cgp_contact_manifdir, paste0(project_id,"_cgp_contact_master_manifest_ids_ony_samp_type.tsv")), header = TRUE, sep="\t", stringsAsFactors = FALSE)
dim(metadata_tab)


#   0.1.3 Merge the samples details and QC table from canapps ------------------
# There is a difference in num of Unique samples and lanes run  
#Then we merge the CGP tables removing the cgp_sp_details Seq.X column 
### There is NO SEQSCAPE manifest so to process the data we will require the following 
cp_qc_details$SANGER_SAMPLE_ID<-lapply(strsplit(cp_qc_details$Sample, split =" ", fixed = T ), function(x){ gsub(")", "", gsub("(", "", x[[2]], fixed=T ), fixed = T) } ) %>% unlist

ttab<-cp_sp_details[match(cp_qc_details$sample, cp_sp_details$Sample), !(colnames(cp_sp_details) %in% c("Seq_X"))   ]
#The merge removes the sample_name from the sample name  details(col1) from samp details table as this is already contained wihtin the QC table
cp_master_tab<-cbind(cp_qc_details , ttab[,c(2:dim(ttab)[2])] ) 
rm(ttab)

#   0.1.4 Load the Sequencescape information  ------------------ 
sqp_manif<-read.csv(file = file.path(seqscape_manifdir, paste0(project_id,"_master_seqscape_manifest.tsv")), header = TRUE, sep="\t", stringsAsFactors = FALSE)
dim(sqp_manif)



#   0.1.4 Merge all the information CGP ans SEqscape into a single manifest ------------------
#get the information form the sequencescape manifest to add into  
ttab<-sqp_manif[ match(cp_master_tab$sample, sqp_manif$SUPPLIER_SAMPLE_NAME) ,]
merged_cp_sqp_manif<- cbind(cp_master_tab, ttab[,c("SANGER_PLATE_ID","WELL","SANGER_SAMPLE_ID","SUPPLIER_SAMPLE_NAME","COHORT","GENDER")])

#   0.1.5 Merge all the metadata and information about what samples are PDXs,  tumours or normal samples  ------------------
merged_cp_sqp_manif <- cbind(merged_cp_sqp_manif, metadata_tab[ match(merged_cp_sqp_manif$sample, metadata_tab$PD_ID), 1:5])


# 0.2.Merge all the information into a manifest ------------------
print("Merge the iRODS manifest with the merged metadata manifest to include the experimental information in a single file ")
merged_cp_sqp_manif<- merged_cp_sqp_manif[match(merged_cp_sqp_manif$sample, raw_manif$sample_supplier_name),]
master_manif<- cbind(raw_manif, merged_cp_sqp_manif[, !(colnames(merged_cp_sqp_manif) %in% c("Project", "Sample", "Comments", "Confirmed.by", "Finalised.by", "sample", "Infuse", "Sample_Issues")) ])
colnames(master_manif)[1]<-"Sanger_sample"


#Print the last manifest
write.table(master_manif,file=file.path(manif_dir,paste(project_id,study_id,"master_manifest_wmetadata_for_proj.txt", sep="_")), quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

