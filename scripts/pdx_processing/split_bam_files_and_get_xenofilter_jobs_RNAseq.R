#!/usr/bin/env Rscript

#/software/team113/dermatlas/R/R-4.2.2/bin/Rscript
###########################################################
#  split_bam_files_and_get_xenofilter_jobs_RNAseq.R
# This script will do the following:
# 


##############################################################
progname<-"split_bam_files_and_get_xenofilter_jobs_RNAseq.R"

#Loading Libraries
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))


#Get the options for tthe script
####### Set the options list of the program
option_list<- list(
  make_option(c("--study_id"), action="store_true", default=NA, type="character", help="This is the sequencescape project ID number "), 
  make_option(c("--manifest"), action="store_true", type="character", default=NA,  help="This is the manifest table generated with fullpath to the manifest table result of PDX_bwa_mem_mapping_jobs_from_master_manif.R "), 
  make_option(c("--projectdir"), action="store_true", type="character", default=NA,  help="This is the fullpath to project directory "), 
  make_option(c("--xfilter_outdir" ), action="store_true", type="character", default=NA, help="This is the location where the final files will be created"),
  make_option(c("--split_nreads" ), action="store_true", type="integer", default=50000000, help="This is the totatl number of reads per splitted file 5e+07(DEFAULT)"),
  make_option(c("--nthreads" ), action="store_true", type="integer", default=4, help="This is number of cores for multithreading tasks")
  
)
#### Indicate that the file follows to the options used
parser<- OptionParser(usage = "%prog [options] --study_id  \"6509\"  --input_tsv  /My/Sequencing_study/project/dir/6509_cohort_indel.marked.target.pass.vep.tsv.gz --outdir /location/to/my/`results/", option_list=option_list)
arguments<- parse_args(parser, positional_arguments = 0) # no files after the options should be provided

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt<-arguments$options
studyid<-opt$study_id
manifest<-opt$manifest
projectdir<-opt$projectdir
xfilter_outdir<-opt$xfilter_outdir
split_nreads<-opt$split_nreads
ncores<-opt$nthreads


#Write some output starting message
write(paste("#########################################\n", progname, date(), "\nThe seqscape_proj_id used was:", studyid,"\n",
            "The manifest provided with the mouse and BAM file information was:",  manifest, "\n",
            "The projectdir provided :",  projectdir, "\n",
            "The xfitler outdir provided:", xfilter_outdir, "\n",
            "The number of  reads to split :",  split_nreads, "\n",
            sep=" "), stdout())

#Verify that the outidr was set
if(is.na(opt$xfilter_outdir)){
  stop("No output directory provided ")
} else {
  if(!dir.exists(xfilter_outdir)){ # If the path doesn't exist create it
    dir.create(xfilter_outdir, recursive = T)
  }
}
#Verify that the study was set 
if(is.na(opt$study_id)){
  stop("No sequencescape Study ID provided  ")
}

#Verify that the manifest file exists
if(is.na(opt$manifest)){
  stop("No manifest  file provided  ")
} else{
  if(!file.exists(opt$manifest)) {
    stop("The manifest provided doesn't exists")
  }
}

#Verify that the outidr was set
if(is.na(opt$xfilter_outdir)){
  stop("No output directory provided ")
} else {
  if(!dir.exists(xfilter_outdir)){ # If the path doesn't exist create it
    dir.create(xfilter_outdir, recursive = T)
  }
}

#Check if samtools is visible 
samtools<-system("samtools --version | head -n 2", intern = T,  wait = TRUE )
if(length(samtools)!=0) {
  write(paste( date(), "\nThe samtools version found was: ", samtools,"\n",
         sep=" "), stdout())
} else {
  stop("Couldn't find samtools ")
}


##############
# STEP1 get filles with all the read names 

#Read the table 
manif<-as.data.frame(data.table::fread(input = manifest,header=TRUE, stringsAsFactors = FALSE))
#Set folders
noddir<- file.path(projectdir, "bams", "NOD_mapping", "NOD_PDXV1")
unfilt_bamdir<-file.path(projectdir, "bams", "RNASEQ_UNFILT")
xfilt_bamdir<-file.path(projectdir, "bams", "RNASEQ_xfilt", "NOD_PDXV1")
scriptsdir<-file.path(projectdir, "scripts", "pdx_processing")
logdir<-file.path(projectdir, "logs")
rsplit_logdir<- file.path(logdir, "read_bam_split")
xfilt_logdir<- file.path(logdir, "xfilt_logdir")
tmpdir<-file.path(projectdir, "tmp")
dir.create(tmpdir, recursive = T)
dir.create(rsplit_logdir, recursive = T)
dir.create(xfilt_logdir, recursive = T)


# First get the read names for all the files
i<-1
new_manif<-NULL
read_name_file<-NULL
hcmds<-NULL
nod_cmds<-NULL
for(i in 1:dim(manif)[1]) {
  temp_bam_basename<-NULL
  temp_rnamfile<-NULL
  temp_sampname<-NULL
  cmd<-NULL
  temp_outdir<-NULL
  temp_hbam<- NULL
  tmp_tot_rfiles<-NULL
  tmp_manif<-NULL
  # Check if the sample needs to be xenofiltered or not
  if(manif$Proc_as_PDX[i]=="N"){
    # Add the sections to the new manifest
      tmp_manif<- cbind(manif[i, ],`tmp_rnam_partfile_list[j]`=NA, tp_hbampartname=NA ) 
      #Add to the new manifest
      tmp_manif$nod_part_bams<- NA
      new_manif<-rbind(new_manif,tmp_manif )
      next
  }
  
 # Get the file names for BAM, output read names and sample being used
  temp_hbam<-manif$unfilt_psample_bam_path[i]
  temp_bam_basename<-basename(manif$unfilt_psample_bam_path[i])
  temp_sampname<-manif$sample[i]
  temp_rnamfile<-file.path(unfilt_bamdir, temp_sampname, gsub(".bam", "_rnames.txt", temp_bam_basename ) )
  temp_rnamfile_basenam<-basename(temp_rnamfile)
  ########### Obtain the read names on the sample bam file
  write( paste0(date(), " Getting the read names for : ", temp_sampname, "\n") , stdout()) 
  #cmd<- paste0("samtools sort --tmpdir ", tmpdir, " -@ ", ncores," -n ", temp_hbam ," | samtools view - | awk '{print $1}' >", temp_rnamfile)
  cmd<- paste0("samtools sort -@ ", ncores," -n ", temp_hbam ," | samtools view - | awk '{print $1}' >", temp_rnamfile)
  write( cmd , stdout()) 
  #system(cmd, wait = T)  
  
  ####### Then split the file, get total number of files to section the read names
  if(manif$total_reads[i]< split_nreads){
    tmp_tot_rfiles<-1 # Total number of read files
  } else {
    tmp_tot_rfiles<-trunc(as.numeric(manif$total_reads[i])/split_nreads) # Total number of read files
  }
  temp_rnamfile_bnam_part<-file.path(unfilt_bamdir, temp_sampname, gsub(".txt","_part_",  temp_rnamfile_basenam, fixed = T) )
  write( paste0(date(), " Getting splittin the read names in  ",tmp_tot_rfiles, " for sample: ",  temp_sampname, "\n") , stdout()) 
  cmd<- paste0("split -n ", tmp_tot_rfiles, " ", temp_rnamfile,  " ", temp_rnamfile_bnam_part, " -a 2 -d ")
  bash_cmd <- paste("bash -c", shQuote(cmd))
  # Print the command to verify
  write(bash_cmd, stdout())
  #system( bash_cmd, wait = T )  
  
  
  #List of partial files 
  tmp_rnam_partfile_list<- system(paste0("ls -1 ", temp_rnamfile_bnam_part ,"*"  ), intern = T)
  ########## Generate the job running to split the BAM files human BAM files
  tmp_hpartbam<-NULL #Temporary  human bam file part vectors
  tmp_nodpartbam<-NULL #Temporary  mouse bam file part vectors
  #Human  BAM splitting jobs
  for( j in 1:length(tmp_rnam_partfile_list)){
    number_part<-NULL
    #Extract the part number
    number_part<-gsub(basename(temp_rnamfile_bnam_part), "", basename(tmp_rnam_partfile_list[j]), fixed = T) # Get the numbers after temp_rnamfile_bnam_part basene
    tp_hbampartname<-NULL #Temporary human bam part name
    tp_hbampartname<- file.path(unfilt_bamdir, temp_sampname, paste0(gsub(".bam", paste0("_part_"),  basename(manif$unfilt_psample_bam_path[i])), number_part, ".bam") )
    # Add the sections to the new manifest
    tmp_manif<- rbind(tmp_manif, cbind(manif[i, ],tmp_rnam_partfile_list[j], tp_hbampartname ) )
    
        #Command for samtools to split the human gene 
    cmd<-NULL
    #cmd<- paste0("samtools view -h -b -tmpdir ", tmpdir, " -@ ", ncores, " --qname-file " , tmp_rnam_partfile_list[j], " ",
    #             manif$unfilt_psample_bam_path[i], " >", tp_hbampartname, " && samtools index ", tp_hbampartname )
    cmd<- paste0("samtools view -h -b -@ ", ncores, " --qname-file " , tmp_rnam_partfile_list[j], " ",
                 manif$unfilt_psample_bam_path[i], " >", tp_hbampartname, " && samtools index ", tp_hbampartname )
    write( paste0(date(), " Splittin the human bam file: ",temp_bam_basename, " for sample: ",  temp_sampname," part :", number_part,  " \n") , stdout()) 
    write(cmd, stdout())
    #Assing job output names 
    sdoutf<-file.path(rsplit_logdir, paste("hum_bam_split_",temp_sampname,"_part", number_part,".o", sep = ""))
    serrf<-file.path(rsplit_logdir, paste("hum_bam_split_",temp_sampname,"_part", number_part,".e", sep = ""))
    #Command for Nodv1 per sample file and index Min 64GBs RAM
    cmd<-as.character(paste("bsub -q normal -M 8000 -R",shQuote("select[mem>8000] rusage[mem=8000] span[hosts=1]") , " -n ", ncores, " -o ", sdoutf," -e ",serrf,
                            " ' ", cmd,
                            " '",
                            sep=""))  
    #Add the command 
    hcmds<- c(hcmds, cmd)
  }
  
  nod_part_bams<-NULL
  ### For NOD BAMs
  for( j in 1:length(tmp_rnam_partfile_list)){
    #NOD  jobs
    number_part<-NULL
    #Extract the part number
    number_part<-gsub(basename(temp_rnamfile_bnam_part), "", basename(tmp_rnam_partfile_list[j]), fixed = T) # Get the numbers after temp_rnamfile_bnam_part basene
    tp_nodbampartname<-NULL #Temporary human bam part name
    tp_nodbampartname<- file.path(noddir, temp_sampname, paste0(gsub(".bam", paste0("_part_"),  basename(manif$bam_path_nodv1[i])), number_part, ".bam") )
    # Add the sections to the new manifest
    nod_part_bams<- c(nod_part_bams,tp_nodbampartname ) 
    
    #Command for samtools to split the NOD BAM
    cmd<-NULL
    #cmd<- paste0("samtools view -h -b -tmpdir ", tmpdir, " -@ ", ncores, " --qname-file " , tmp_rnam_partfile_list[j], " ",
    #             manif$bam_path_nodv1[i], " >", tp_nodbampartname, " && samtools index ", tp_nodbampartname )
    cmd<- paste0("samtools view -h -b -@ ", ncores, " --qname-file " , tmp_rnam_partfile_list[j], " ",
                 manif$bam_path_nodv1[i], " >", tp_nodbampartname, " && samtools index ", tp_nodbampartname )
    
    write( paste0(date(), " Splittin the NOD bam file: ", basename(manif$bam_path_nodv1[i]), " for sample: ",  temp_sampname," part :", number_part,  " \n") , stdout()) 
    write(bash_cmd, stdout())
    #Assing job output names 
    sdoutf<-file.path(rsplit_logdir, paste("nod_bam_split_",temp_sampname,"_part", number_part,".o", sep = ""))
    serrf<-file.path(rsplit_logdir, paste("nod_bam_split_",temp_sampname,"_part", number_part,".e", sep = ""))
    #Command for Nodv1 per sample file and index Min 64GBs RAM
    cmd<-as.character(paste("bsub -q normal -M 8000 -R",shQuote("select[mem>8000] rusage[mem=8000] span[hosts=1]") , " -n ", ncores, " -o ", sdoutf," -e ",serrf,
                            " ' ", cmd,
                            " '",
                            sep=""))  
    
    #Add the command 
    nod_cmds<- c(nod_cmds, cmd)
  }
  
  #Add to the new manifest
  tmp_manif$nod_part_bams<- nod_part_bams
  
  new_manif<-rbind(new_manif,tmp_manif)
}


########## Write the sh file that will submit the jobs to split the BAM files
#Generate the .sh files with the submissions for the BAM splitting jobs of human  
write.table(c("#!/bin/bash", hcmds, nod_cmds ), file=file.path(scriptsdir, paste(studyid , "_bamsplittin_by_read_names_jobs.sh", sep="")), quote = F, col.names = F, row.names = F, sep = "\n")


#Fixed Column names for unfilter part bams
colnames(new_manif)<-gsub("tmp_rnam_partfile_list[j]", "read_names", colnames(new_manif), fixed = T)
colnames(new_manif)<-gsub("tp_hbampartname", "hum_unfilt_partbam", colnames(new_manif), fixed = T)



### Generate the jobs for Xenofilter for the BAMs
write( paste0(date(), " Generating the Xenofilter job commands for all the samples  \n") , stdout()) 
xfilt_bam_psample_path_nodv1_part<-NULL
xfilt_part_cmds<-NULL
xenofilter_script<-file.path(scriptsdir, "bam_Xenofilter_rg.R")
i<-1
for (i in 1:dim(new_manif)[1] ){
  # Get the commands for submitting Xenofilter Jobs
  temp_sampname<-NULL
  cmd<-NULL
  temp_outdir<-NULL
  tmp_manif<-NULL
  number_part<-NULL
  temp_part_xfiltbamdir<-NULL
  
  # Check if the sample needs to be xenofiltered or not
  if(new_manif$Proc_as_PDX[i]=="N"){
    # Add the value for this instances
    temp_part_xfiltbamdir<-NA
    xfilt_bam_psample_path_nodv1_part<-c(xfilt_bam_psample_path_nodv1_part, temp_part_xfiltbamdir )
    next
  }
  
  #Variables
  temp_sampname<-new_manif$sample[i]
  #Extract the part number
  number_part<-unlist(strsplit(basename(new_manif$read_names[i]), "_part_", fixed = T))[2]
  #Output name of the part X filt file
  temp_part_xfiltbamdir<-file.path(xfilter_outdir, paste0(temp_sampname, "_",  number_part),"Filtered_bams",
                                   paste0(temp_sampname,"_",number_part, ".sample.dupmarked.mXfilt_Filtered.bam") )
  #Command for submitting the Xenofilter job commands 
  cmd<-NULL
  cmd<- paste0("Rscript ", xenofilter_script, 
               " --sample_name ", temp_sampname, 
               " --human_bam ", new_manif$hum_unfilt_partbam[i],
               " --mouse_bam " , new_manif$nod_part_bams[i],
               " --rg_id ", number_part, 
               " --outdir " , xfilter_outdir,
               " --ncpu 1"  
                )
  #Assing job output names 
  sdoutf<-file.path(xfilt_logdir, paste("xenofilter_",temp_sampname,"_part", number_part,".o", sep = ""))
  serrf<-file.path(xfilt_logdir, paste("xenofilter_",temp_sampname,"_part", number_part,".e", sep = ""))
  #Command for Xenofilter per sample file and index Min 36GBs RAM for 50e+6 Reads
  cmd<-as.character(paste("bsub -q long -M 36000 -R",shQuote("select[mem>36000] rusage[mem=36000] span[hosts=1]") , " -n ", ncores, " -o ", sdoutf," -e ",serrf,
                          " ' ", cmd,
                          " '",
                          sep=""))  
  #xfilter commands
  xfilt_part_cmds<- c(xfilt_part_cmds, cmd)
  # Name of the files
  xfilt_bam_psample_path_nodv1_part<-c(xfilt_bam_psample_path_nodv1_part, temp_part_xfiltbamdir )
}

new_manif$xfilt_bam_psample_path_nodv1_part<- xfilt_bam_psample_path_nodv1_part

#Generate the .sh files with the submissions for the BAM splitting jobs of human  
write.table(c("#!/bin/bash", xfilt_part_cmds), file=file.path(scriptsdir, paste(studyid , "_Xenofilter_by_read_names_parts_jobs.sh", sep="")), quote = F, col.names = F, row.names = F, sep = "\n")

############ Generate the jobs to merge the BAM files
write( paste0(date(), " Generating the job commands  to merge Xenofiltered files samples \n") , stdout()) 
xfilt_merge_cmds<-NULL
xfilt_bam_psample_path_nodv1<-NULL
sample_list<- unique(new_manif$sample)
for (i in 1:length(sample_list) ){
  # Get the commands for submitting Xenofilter Jobs
  temp_sampname<-NULL
  temp_smanif<-NULL
  cmd<-NULL
  temp_merged_xfile<- NULL
   
  # Check if the sample needs to be xenofiltered or not
  if(new_manif$Proc_as_PDX[i]=="N"){
    # Add the value for this instances
    temp_merged_xfile<-NA
    xfilt_bam_psample_path_nodv1<-c(xfilt_bam_psample_path_nodv1, temp_merged_xfile )
    next
  }
  
  #Variables
  temp_sampname<-sample_list[i]
  temp_smanif<- new_manif[new_manif$sample == temp_sampname, ]
  #Output name of the Xfilt file
  dir.create(file.path(xfilter_outdir, paste0(temp_sampname),"Filtered_bams"), recursive=T)
  temp_merged_xfile<-file.path(xfilter_outdir, paste0(temp_sampname),"Filtered_bams", 
                               paste0(temp_sampname,".sample.mXfilt_Filtered.bam") )
    #Command for submitting the Xenofilter job commands 
  cmd<-NULL
  cmd<- paste0( paste("samtools merge --threads", ncores, 
              "-o", temp_merged_xfile,
              paste(temp_smanif$xfilt_bam_psample_path_nodv1_part, collapse = " "), sep = " " ), " ; ", 
              paste0("samtools index ", temp_merged_xfile))
              
  #Assing job output names 
  sdoutf<-file.path(xfilt_logdir, paste("xenofilter_merge_",temp_sampname,".o", sep = ""))
  serrf<-file.path(xfilt_logdir, paste("xenofilter_merge_",temp_sampname,".e", sep = ""))
  #Command for Xenofilter per sample file and index Min 36GBs RAM for 50e+6 Reads
  cmd<-as.character(paste("bsub -q normal -M 8000 -R",shQuote("select[mem>8000] rusage[mem=8000] span[hosts=1]") , " -n ", ncores, " -o ", sdoutf," -e ",serrf,
                          " ' ", cmd,
                          " '",
                          sep=""))  
  #xfilter commands
  xfilt_merge_cmds<- c(xfilt_merge_cmds, cmd)
  # Name of the files
  xfilt_bam_psample_path_nodv1<-c(xfilt_bam_psample_path_nodv1, rep(temp_merged_xfile, dim(temp_smanif)[1] ))
}

new_manif$xfilt_bam_psample_path_nodv1<- xfilt_bam_psample_path_nodv1

#Generate the .sh files with the submissions for the BAM splitting jobs of human  
write.table(c("#!/bin/bash", xfilt_merge_cmds), file=file.path(scriptsdir, paste(studyid , "_Xenofilter_merged_filtered_parts_jobs.sh", sep="")), quote = F, col.names = F, row.names = F, sep = "\n")

write.table( new_manif, file.path(projectdir, "metadata", "manifests", paste0(studyid, "_cram_manifest_INFO_from_iRODS_all_wbam_counts_qc_wNOD_bams_part.txt")), quote = F, col.names = T, row.names = F, sep = "\t")


