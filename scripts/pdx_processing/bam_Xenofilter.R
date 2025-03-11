#!/software/R-4.1.0/bin/Rscript 
#################################################################
# bam_Xenofilter.R
# This script takes a sample_name, human mapped bamfile, mouse mapped bamfile, output directory, number of CPUS
# to run Xenofilter and filter the mouse reads. 
# This script will create the following:
#   a) Create a directory within the provided output directorey  with the sample_name provided
#   b) Filter the provided bam file with XenofilteR 
#   c) Create a tab separated file with the following information c("sample", "unfilt_bamfile", "filtered_read_pairs", "total_read_pairs", "percent_filtered")
# 
# Created: 08/2022  
#Author: Martin Del Castillo Velasco-Herrera - mdc1@sanger.ac.uk
###################################################################
progname<-"bam_Xenofilter.R"
print("Setting up the environment ")
# Check for requirements

#Loading Libraries
suppressMessages(library(XenofilteR))
suppressMessages(library(pryr))
suppressMessages(library(optparse))

#Get the options for tthe script
####### Set the options list of the program
option_list<- list(
  make_option(c("--sample_name"), action="store_true", default=NA, type="character", help="This is the name of the sample being analysed and Prefix for the ouput file"), 
  make_option(c("--human_bam" ), action="store_true",type="character", default=NA, help="Full path to the BAM file with reads mapped to the human genome"), 
  make_option(c("--mouse_bam" ), action="store_true",type="character", default=NA, help="Full path to the BAM file with reads mapped to the human genome"),
  make_option(c("--outdir" ), action="store_true",type="character", default=NA, help="Full path to the directory to write the ouput files"),
  make_option(c("--ncpu" ), action="store_true",type="character", default="1", help="Number of CPUs (workers) used for SnowParams")
)
#### Indicate that the file follows to the options used
parser<- OptionParser(usage = paste(progname, " [options] --sample_name PDXXXXXa --human_bam  /PATH/TO/PDXXXXa.human.bam --mouse_bam  /PATH/TO/PDXXXXa.mouse.bam --outdir /PATH/TO/OUTPUT --ncpu 1 ",sep=""), option_list=option_list)
arguments<- parse_args(parser, positional_arguments = 0) # no files after the options shoudl be provided
#FOR testing purposes ONLY
# arguments<- parse_args(parser,
#                        args=c("--sample_name","PD53330a",
#                               "--human_bam","/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/BAMS/WES/PD53330a.sample.dupmarked.bam",
#                               "--mouse_bam","/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/BAMS/NOD_mapping/NOD_PDXV1/PD53330a/PD53330a_NOD_PDXV1.aln.sort.out.bam",
#                               "--outdir", "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/test_xenofilter",
#                               "--ncpu", "1"
#                               ), #For tests
#                        # args=c("--manifest","6591_2775_master_manifest_wmetadata_for_proj.txt","--projectdir","/Users/mdc1/Desktop/team113sc124/projects/6591_PDX_models_Latin_America_RNAseq"), #For tests from PC
#                        positional_arguments = 0) # no files after the options should be provided

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt<-arguments$options

# Check the PARAMETERS provided ----- 
#Verify that the outidr was set and create it if it doesn't exists
if(is.na(opt$outdir)){
  outdir<-getwd()
} else {
  outdir<-file.path(opt$outdir, opt$sample_name)
}
if(!dir.exists(outdir)){ # Create the outdir if it doesn't exist
  dir.create(outdir, recursive = T )
}
#Exit if the sample name is not provided
if(is.na(opt$sample_name)){
  stop(paste("NO sample_name provided, provide a sample_name -required_parameter-", sep=""))
}
#Verify that the human bam is provided and exists
if(is.na(opt$human_bam)){
  stop(paste("NO human_bam file provided, provide a human_bam  -required_parameter-", sep=""))
} else if(!file.exists(opt$human_bam)){
  stop(paste("The human_bam provided: ", opt$human_bam, " doesn't exists!", sep=""))
}
#Verify that the mouse bam is provided and exists
if(is.na(opt$mouse_bam)){
  stop(paste("NO mouse_bam file provided, provide a mouse_bam  -required_parameter-", sep=""))
} else if(!file.exists(opt$mouse_bam)){
  stop(paste("The mouse_bam provided: ", opt$mouse_bam, " doesn't exists!", sep=""))
}

#Write some output starting message
write(paste("#########################################\n", progname, date(),"\nThe sample to be proccessed is: ",opt$sample_name, "\nThe human_bam provided: ", opt$human_bam,"\n",
            "The mouse_bam provided: ", opt$mouse_bam, "\n", "The output dir is: ", opt$outdir, "\n",
            "The samples will be filtered using nCPUs:", opt$ncpu, sep=""), stdout())

#######################################################################
#Defining the directories for the analyses and GENERAL VARIABLES ------

#Set of name of the variables to use
sample_name<-opt$sample_name
hbam<-opt$human_bam
mbam<-opt$mouse_bam

#Create the sample.list data frame
outfname<-paste0(sample_name, ".sample.dupmarked.mXfilt.bam")
sample.list<-data.frame( hbam, mbam )

#Define the environment variables for the BiocParallel
# SnowParam : used for computing in distributed memory
# MulticoreParam : used for computing in shared memory
# BatchtoolsParam: for computing with cluster schedulers
# DoparParam: for computing with foreach
# bp.param<- SnowParam(workers=opt$ncpu, type="SOCK")
# bp.param<- MulticoreParam(workers=opt$ncpu)
bp.param<- SnowParam(workers=as.numeric(opt$ncpu), type="SOCK")


# Run Xenofilter -------
write(paste("\n---------------------------------------\n", progname, date(),"\n Running Xenofilter for : ",opt$sample_name, "\n", sep=""), stdout())
print(mem_used())
system.time(
  tvar<-capture.output(XenofilteR(sample.list, destination.folder = file.path(outdir) , bp.param = bp.param, output.names=outfname, MM_threshold=4 ), stdout())
)



# Xenofilter Filter stats -------
#Get the statistics of filtered numbers 
tstats<-tvar[grepl("  - Filtered ", tvar, fixed = T)]
#Remove the additional symbols from the text whose output is something like :
# [1] \"Test_hg19_NRAS.bam  - Filtered 59 read pairs out of 645  -  9.15 Percent\""
tstats<- gsub("[1] \"", "", tstats, fixed = T) #remove the beginning 
tstats<- gsub("  - Filtered ", "\t", tstats, fixed = T) #replace   - Filtered  to \t
tstats<- gsub(" read pairs out of ", "\t", tstats, fixed = T) #replace  "read pairs out of "  for \t
tstats<- gsub("  -  ", "\t", tstats, fixed = T) #replace  "  -  "  for \t
tstats<- gsub(" Percent\"", "", tstats, fixed = T) #replace  " Percent\""  for
#Removing the tabs to get a vector and then add the sample name to create a vector 
tstats<-strsplit(tstats, split = "\t", fixed = T)[[1]]
tstats<- c(sample_name, tstats)
names(tstats)<- c("sample", "unfilt_bamfile", "filtered_read_pairs", "total_read_pairs", "percent_filtered")

#Save the stats in the same folder 
write.table(data.frame(t(tstats)), file=file.path(outdir, "Filtered_bams", gsub(".bam", "_Filtered.bam.xfiltstats", outfname, fixed = T)), row.name=FALSE, col.names=TRUE, quote=FALSE,  sep="\t")  


