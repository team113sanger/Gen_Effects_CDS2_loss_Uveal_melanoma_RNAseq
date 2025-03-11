#!/bin/bash

#Clean the environment
module purge 

# R and R libraries 
_R_BIN_PATH="/software/team113/dermatlas/R/R-4.2.2/bin"
export PATH="${_R_BIN_PATH:?empty-path-variable}:${PATH}"
export R_LIBS="/lustre/scratch124/casm/team113/projects/7688_7687_Gen_Effects_CDS2_loss_Uveal_melanoma/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq/renv/library/R-4.2/x86_64-pc-linux-gnu"


#load iRODS module 
module load IRODS/1.0

# Required samtools versions 
module load samtools-1.14/python-3.12.0
module load star/2.7.10a 

NOD_PDXV1_REFDIR="/lustre/scratch124/casm/team113/projects/7688_7687_Gen_Effects_CDS2_loss_Uveal_melanoma/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq/reference/NOD_ShiLtJ_V1_PDX_ref"
HUMSQE_REF="/lustre/scratch125/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa"


##HTseq-counts version 0.13.5
#"/software/team113/users/mdc1/python/bin/htseq-count"
