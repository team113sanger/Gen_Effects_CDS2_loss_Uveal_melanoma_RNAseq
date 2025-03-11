#!/usr/bin/env Rscript
###########################################################
#  Build_manifest_from_irods_cram_information.R
# This script will do the following:
#  1.Get the list of cram files for the specified sequencescape project from iRODS
#  2.Get the metadata for each CRAM file from iRODs to create a sample manifest
#
# Author:Martin Del Castillo Velasco Herrera
##############################################################
progname<-"Build_manifest_from_irods_cram_information.R"
#Loading Librarie
suppressMessages(library(optparse))

#Get the options for tthe script
####### Set the options list of the program
option_list<- list(
		     make_option(c("--seqscape_proj_id"), action="store_true", default=NA, type="character", help="This is the sequencescape project ID number "), 
		       make_option(c("--outdir" ), action="store_true",type="character", default=NA, help="This is the location where the file will be created(if provided else the current directory)")
		     )
#### Indicate that the file follows to the options used
parser<- OptionParser(usage = "%prog [options] --seqscape_proj_id  \"My project title\"", option_list=option_list)
arguments<- parse_args(parser, positional_arguments = 0) # no files after the options shoudl be provided

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt<-arguments$options
study<-opt$seqscape_proj_id

#Verify that the outidr was set
if(is.na(opt$outdir)){
	  outdir<-paste(getwd(),"/",sep="")
} else {
	  outdir<-paste(opt$outdir,"/",sep="")
}
#Verify that the study was set 
if(is.na(opt$seqscape_proj_id)){
	  stop("No sequencescape Study ID provided  ")
}

#Write some output starting message
write(paste("#########################################\n", progname, date(), "\nThe seqscape_proj_id used was:", study,"\n",
	                "The outdir with the manifest is the following:", outdir, "\n", sep=" "), stdout())

#############################################################################################
#  DEFINE the functions of the script
#############################################################################################
#To get the manifest from a project :
#To get the information from a cram file in iRODS run:
#imeta ls -d /seq/21929/21929_1#223.cram | grep
#The important CRAM file metadata is the following:
# sample - SANGER_sample_ID
# 4458STDY6706820
# sample_supplier_name - supplier sample _name
# SCGC--1342_E04
# sample_donor_id  - SANGER_sample_ID
# 4458STDY6706820
# sample_accession_number - EGA - ENA ACCESSION NUMBER
# EGAN00001497313
# study_accession_number - ENA STUDY accession number
# EGAS00001002104
# sample_common_name - Species
# Homo sapiens
# md5  
# a440bb5848b61d5789436b3f651ca2ddid_run
# id_run -id run
# 21929
# library_id - library id
# 18474284
# total_reads - total reads generated
# 1562728
# study_id - Sequencescape study ID
# 4458
# study_title - Sequencescape study title
# value: Psoriatic arthritis
# reference - reference genome used
# value: /lustre/scratch117/core/sciops_repository/references/Homo_sapiens/1000Genomes_hs37d5/all/bwa0_6/hs37d5.fa
# library_type - library title
# value: Nextera Dual Index qPCR only
###############################
# get_cram_info_item_value_from_IRODs - FUNCTION that will extract the metadata value desired form irods for a cram file
###############################
get_cram_info_item_value_from_IRODs<- function(cram,value){
	  if((!is.character(value)) & (!is.vector(value)) ){
		      stop("The value is not a character vector")
  }
  if( (!is.character(cram)) & (length(value)!=1) ){
	      stop("cram is not a character vector of length 1")
    }
    #VARIABLEs for the function
    temp_imeta<-NULL
    temp_metadata<-NULL
      #Make a system call to get the information
      # imeta ls -d /seq/21929/21929_1#223.cram | grep
      # cram<-"/seq/21929/21929_1#223.cram"
      #values<-c("sample","sample_supplier_name ","sample_donor_id","sample_accession_number","study_accession_number","sample_common_name","md5  ","id_run","library_id","total_reads","study_id","study_title","reference","library_type")
      # GEt the cram info
      temp_imeta<- system( paste("imeta ls -d ", cram, sep=""), intern = T)
      #extract the values 
      for (i in 1:length(value)){
	          #Get the position of the matching attribute information
	          pos<-match(paste("attribute: ",value[i], sep=""),temp_imeta)
          # Get the value of the attribute
          temp_metadata<-c(temp_metadata,gsub("value: ", "", temp_imeta[pos+1]))
	    }
        names(temp_metadata)<- value
        #
        return(temp_metadata)
}

########################################
#  Function to get the cram file list from iRODS
########################################
get_cram_list_from_irods<-function(study_id){
	  if( (!is.character(study_id)) & (length(study_id)!=1) ){
		      stop("cram is not a character vector of length 1")
  }
  #NEcessary variables
  temp_collection<-NULL
    temp_object<- NULL
    temp_cram_list<-NULL
      #Get the study collection of folders from imeta e.g imeta qu -z seq -d study = "Psoriatic arthritis" and type = "cram" and target = 1 | grep “collection:" >collection
      # study_id<-"Psoriatic arthritis"
      print(paste("Retrieving Project ID: ",study_id, " using :", sep=""))
      print(paste("imeta qu -z seq -d study_id = \"",study_id,"\" and type = \"cram\" and target = 1| grep collection", sep="" ))
        temp_collection<-gsub("collection: ", "",system(paste("imeta qu -z seq -d study_id = \"",study_id,"\" and type = \"cram\" and target = 1| grep collection", sep="" ),intern = T))
        temp_object<-gsub("dataObj: ", "",system(paste("imeta qu -z seq -d study_id = \"",study_id,"\" and type = \"cram\" and target = 1 | grep dataObj", sep="" ),intern = T))
	  #Get_cramlist
	  temp_cram_list<- paste(temp_collection, temp_object,sep="/")
	  return(temp_cram_list)
}

########################################
#  Function to get the cram file names from the cramfile_list
########################################
get_cram_fname_from_cramlist<-function(crams){
	  if((!is.character(crams)) & (!is.vector(crams)) ){
		      stop("The value is not a character vector")
  }
  #Variables for the fucntion
  cram_name<-NULL
    for(i in 1:length(crams)){
	        temp<-NULL
      temp<-tail(unlist(strsplit(crams[i], "/")), n=1)
          cram_name<-c(cram_name,temp )
        }
    return(cram_name)
}

###########################################################################
# Functions to get the C>A information from the substitution metrics list
###########################################################################
get_CtoA_results<-function(subsf){
  if((!is.character(subsf)) & (!is.vector(subsf)) ){
    stop("The value is not a character vector")
  }
  #Variables for the funrction
  CtoA_vals<-NULL
  for(i in 1:length(subsf)){
    temp<-NULL
    cmd<-NULL
    cmd<-paste("iget ", subsf[i], " - | grep 'art_oxH' | cut -f 2 -d ' ' ")
    temp<-system(command=cmd , intern = TRUE)
    CtoA_vals<-c(CtoA_vals, temp)
  }
  return(CtoA_vals)
}

################################################################
#  1.Get the list of cram files for the project from iRODS
################################################################
#For testing purposes
# study<-"6633"
# outdir<- "/lustre/scratch119/realdata/mdt1/team113/projects/6633_PDX_models_Latin_America_WES/manifests"
# study_name<-"Psoriatic_arthritis"
study_name<-gsub(" ", "_", study)
#info from each cram file to get
values<-c("sample","sample_supplier_name","sample_donor_id","sample_accession_number","study_accession_number","sample_common_name","md5","id_run","lane","is_paired_read","tag_index","library_id","total_reads","study_id","study_title","reference","library_type")
#get cramfile list
cram_list<-get_cram_list_from_irods(study)

################################################################
#  2.Get the metadata for each CRAM file from iRODs to create a sample manifest
################################################################
manifest<-NULL
for(j in 1:length(cram_list)) {
	  temp<-NULL
  temp<-get_cram_info_item_value_from_IRODs(cram_list[j], values)
    manifest<-rbind(manifest,temp)
    if((j %% 100)==0){
	        print(paste("getting info for CRAM number ", j, sep=""))
      }
}  
#Add sensible colnames
colnames(manifest)<-values
row.names(manifest)<- 1:dim(manifest)[1]
#Transform to a data fram
manifest<-as.data.frame(manifest,stringsAsFactors=F)
#Add the cramfile names and irods location
manifest$cram<-get_cram_fname_from_cramlist(cram_list)
manifest$cram_irods_location<- cram_list
#Write the manifest after ordering by 
manifest<-manifest[order(manifest$sample_supplier_name), ]
#This will generate the information for the C>A file name in the form the name of the CRAM file 
manifest$ca_metrics_file<- gsub(".cram", ".substitution_metrics.txt", manifest$cram_irods_location ) #The file name
#Call the function to add the 
manifest$art_oxH<-get_CtoA_results(manifest$ca_metrics_file)

# 1.4 Print the violin plot of insert size  ------------------
manifest$art_oxH<-as.numeric(manifest$art_oxH)
#Violin plote with value for the cohort
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
viol_ctoAplot<- ggplot(reshape2::melt(manifest, id.vars="sample_supplier_name", measure.vars="art_oxH"), aes(x=variable, y=value))+
  geom_violin(aes(fill=variable))+
  scale_fill_manual(values = "lightblue")+ # to fill the violing plot 
  ggtitle(paste(study_name, "\n C>A art_oxH score distribution "))+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.text.y=element_text(colour="Black",size =13),
    axis.title.x=element_text(face="bold", size=14),
    axis.title.y=element_text(face="bold", size=14),
    legend.text = element_text(face="bold", size=11),
    plot.title = element_text(hjust = 0.5))+
  ylab("C>A art_oxH")+
  xlab("Cohort")+
  geom_boxplot(width=0.1)
ggsave(filename = file.path(paste(outdir,study_name,"_CtoA_vplot.pdf",sep="")), 
       viol_ctoAplot)

write.table(manifest, file = paste(outdir,study_name,"_cram_manifest_INFO_from_iRODS.txt",sep=""), quote=F, row.names=FALSE, sep="\t")  


# ################################################################
# #  3. Generate the summary of the sample data 
# ################################################################
# #Get the unique sample names
unique_sample_sup_names<- unique(manifest$sample_supplier_name)
write.table(unique_sample_sup_names, file = paste(outdir,study_name,"_manifest_uniq_sup_sample_names.txt",sep=""), quote=F, col.names=FALSE,row.names=FALSE, sep="\t")  


