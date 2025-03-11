#!/bin/bash
module purge
module load rocker/rver/4.4.1
module load git

# For R 4.4.1
_R_BIN_PATH="/software/isg/containers/R/rver/bin/R"
export PATH="${_R_BIN_PATH:?empty-path-variable}:${PATH}"
export R_LIBS="/lustre/scratch124/casm/team113/projects/7688_7687_Gen_Effects_CDS2_loss_Uveal_melanoma/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq/scripts/expression/renv/library/linux-ubuntu-jammy/R-4.4/x86_64-pc-linux-gnu"


