#!/usr/bin/env Rscript
##########################################################
#  collate_RNA_xfilter_stats_from_manif_and_plots.R
# This script will do the following:
# 

##############################################################
progname<-"collate_RNA_xfilter_stats_from_manif_and_plots.R"

#Loading Libraries
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))


#Get the options for tthe script
####### Set the options list of the program
option_list<- list(
  make_option(c("--study_id"), action="store_true", default=NA, type="character", help="This is the sequencescape project ID number "), 
  make_option(c("--manifest"), action="store_true", type="character", default=NA,  help="This is the manifest table generated with fullpath to the manifest table result of PDX_bwa_mem_mapping_jobs_from_master_manif.R "), 
  make_option(c("--projectdir"), action="store_true", type="character", default=NA,  help="This is the fullpath to project directory "), 
  make_option(c("--xfilter_outdir" ), action="store_true", type="character", default=NA, help="This is the location where the final files will be created")
)
#### Indicate that the file follows to the options used
parser<- OptionParser(usage = "%prog [options] --study_id  \"7688\"  --manifest  /My/Sequencing_study/project/dir/7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt --projectdir /My/Sequencing_study/project/dir --xfilter_outdir /My/Sequencing_study/project/dir/analysis/bam_xfilterstats", option_list=option_list)
arguments<- parse_args(parser, positional_arguments = 0) # no files after the options should be provided

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt<-arguments$options
studyid<-opt$study_id
manifest<-opt$manifest
projectdir<-opt$projectdir
xfilter_outdir<-opt$xfilter_outdir


#Write some output starting message
write(paste("#########################################\n", progname, date(), "\nThe seqscape_proj_id used was:", studyid,"\n",
            "The manifest provided with the mouse and BAM file information was:",  manifest, "\n",
            "The projectdir provided :",  projectdir, "\n",
            "The xfitler outdir provided:", xfilter_outdir, "\n",
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


##############
# STEP1  read the manifest and .

#Read the table 
manif<-as.data.frame(data.table::fread(input = manifest,header=TRUE, stringsAsFactors = FALSE))

# Remove samples that are not marked as processed as PDX
manif<- manif[ manif$Proc_as_PDX=="Y",]

# merge table
xftab<-NULL
for(i in 1:dim(manif)[1]) {
  temp<-NULL
  tempfname<-paste0(manif$xfilt_bam_psample_path_nodv1_part[i], ".xfiltstats")
  temp<-read.csv(file=tempfname, header=TRUE, stringsAsFactors=F, sep="\t")
  xftab<-rbind(xftab,temp )
}

# Write the table with the resulting parameters for all the part samples
xfiltab_part_fname<-paste0(studyid,"_cohort_per_part.xfiltstats"  )
fwrite( xftab, file=file.path(xfilter_outdir, xfiltab_part_fname ), quote=FALSE, row.names = FALSE, col.names=TRUE, sep="\t" )

# STEP2 collate data by sample summary and save in final tab and per folder
samplesid<-unique(manif$sample)

psamp_xfatb<-NULL
for(i in 1:length(samplesid)){
  tmptab<-NULL
  n_parts<-NULL
  tftab<- NULL 
  tmpxstats_file<-NULL
  tmptab<- xftab[xftab$sample==samplesid[i],]
  tftab<- tmptab[1,]
  n_parts<- length(unique(tmptab$unfilt_bamfile))
  # Addign the values fof the sample
  #Get the name of the bam file per sample
  tftab$unfilt_bamfile<- unique(manif$unfilt_psample_bam_path[manif$sample==samplesid[i]])
  tftab$filtered_read_pairs<- sum(tmptab$filtered_read_pairs)
  tftab$total_read_pairs<- sum(tmptab$total_read_pairs)
  tftab$percent_filtered<- ((tftab$filtered_read_pair/tftab$total_read_pairs)*100)
  #Generate the table in both the BAM folder and add it to the final one for results 
  tmpxstats_file<- paste0(unique(manif$xfilt_bam_psample_path_nodv1[manif$sample==samplesid[i]]), ".xfiltstats")
  fwrite(tftab, tmpxstats_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  #Add information to cohort per sample table
  psamp_xfatb<- rbind(psamp_xfatb, tftab)
}

# Write the cohort table
xfiltab_psample_fname<-paste0(studyid,"_cohort_per_sample.xfiltstats"  )
fwrite( psamp_xfatb, file=file.path(xfilter_outdir, xfiltab_psample_fname ), quote=FALSE, row.names = FALSE, col.names=TRUE, sep="\t" )


# Step 3 plot the observed filterstats

vcol<- sqr.col<-RColorBrewer::brewer.pal(n=12, name="Set3")[1]
vpl<- ggplot(reshape2::melt(psamp_xfatb, id.vars="sample", measure.vars="percent_filtered"), aes(x=variable, y=value))+
  geom_violin(aes(fill=variable))+
  # scale_fill_manual(values = "lightblue")+ # to fill the violing plot 
  scale_fill_manual(values = sqr.col)+ 
  ggtitle(paste0("Percentage of mouse reads filtered\n ", studyid))+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.text.y=element_text(colour="Black",size =13),
    axis.title.x=element_text(face="bold", size=14),
    axis.title.y=element_text(face="bold", size=14),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10))+
  ylab("Percentage of mouse reads filtered")+
  xlab("Samples")+
  geom_boxplot(width=0.1)

ggsave(filename= file.path(xfilter_outdir, paste0(studyid, "_per_sample_pdx_xfiltstats_vpl.pdf")), 
       vpl)


bpl<- ggplot(reshape2::melt(psamp_xfatb, id.vars="sample", measure.vars="percent_filtered"), aes(x=sample, y=value))+
  geom_bar( stat = "identity",  position = "dodge", aes(fill=variable))+
  scale_fill_manual(values = sqr.col)+ 
  ggtitle(paste0("Percentage of mouse reads filtered\n ", studyid))+
  theme_bw()+
  theme(
    axis.text.x=element_text(angle = 90 , size = 9,vjust = 0.5, hjust=1 ),
    axis.text.y=element_text(colour="Black",size =13),
    axis.title.x=element_text(face="bold", size=14),
    axis.title.y=element_text(face="bold", size=14),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 10))+
  ylab("Percentage of mouse reads filtered")+
  xlab("Samples")
ggsave(filename= file.path(xfilter_outdir, paste0(studyid, "_per_sample_pdx_xfiltstats_bpl.pdf")), 
       bpl)

