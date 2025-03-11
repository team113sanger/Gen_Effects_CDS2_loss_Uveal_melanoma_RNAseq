# Transcriptional quantification using Kallisto

## Overview

This document describes the steps taken to quantify the abundances of transcripts on the RNAseq data generated from the transcriptome of all the samples. The data from both _CDS2_ targeting experiments was quantified using `kallisto 0.51.1` and the ENSEMBLv103 annotation of the `GRCh38` human reference genome. The quantification was performed using the `kallisto quant` function.

All the scripts and code mentioned below can be found in the `scripts/expression` directory.

## Required Environment variables and software

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.2.2`[**here**](https://cran.r-project.org/)
- **singularity** version `3.11.4` [**here**](https://sylabs.io/singularity/)
- **samtools**: `v1.14` [**here**](https://github.com/samtools/samtools)
- **kallisto:** kallisto `0.51.1` [**here**](https://github.com/pachterlab/kallisto/tree/v0.51.1)
-Dataset
  - ENSEMBL version `103`[**here**](http://feb2021.archive.ensembl.org)

### Generation of the Kallisto singularity image

The instructions of how to generate a singularity image that contains `kallisto 0.51.1` can be found in the [Kallisto_image_obtention.md](../scripts/singularity_images/kallisto/Kallisto_image_obtention.md) file.

### R environment

If you're interested in reproducing the R environment, for a the code used in R 4.4.1 run the following commands within `R v4.4.1`, change the path on `projectdir` to the path where the repository was cloned into:

```R
projectdir<-"/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq"
exp_processing_dir<- file.path(projectdir,"scripts/expression")

setwd(pdx_processing_dir)
install.packages("renv")
library(renv)
renv::restore(lockfile=file.path(exp_processing_dir, "renv.lock")) # To rebuild an environment from the renv.lockfile

```


## Kallisto index generation for GRCh38 ENSEMBL v103 annotation

To generate the index with the sequence of the transcripts in ENSEMBLs v103 annotation of the GRCh38 human reference genome with `kallisto index`, the following steps were taken:

**NOTE**
Set the `PROJECTDIR` variable to the path where the repository was cloned into 

```bash
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq
STUDY=7687
cd ${PROJECTDIR}

# Directory to store the reference files
ENSV103_DIR=${PROJECTDIR}/reference/kallisto/GRCh38_ensv103
mkdir -p ${ENSV103_DIR}
## KAllisto singularity image location
KALLISTO_DIR=${PROJECTDIR}/scripts/singularity_images/kallisto/


# Download the required files Human
#Genome: 
cd ${ENSV103_DIR}
wget https://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ./
#GTF: 
wget https://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz
# CDNA FASTAS
wget https://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# NC RNA FASTAS
wget https://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
zcat ${ENSV103_DIR}/Homo_sapiens.GRCh38.cdna.all.fa.gz ${ENSV103_DIR}/Homo_sapiens.GRCh38.ncrna.fa.gz | gzip - >${ENSV103_DIR}/Homo_sapiens.GRCh38.103.rna.fa.gz

# Run  kallisto index
singularity exec -e -B `pwd` -B / -B ${PROJECTDIR:?unset}  \
${KALLISTO_DIR:?unset}/kallisto_0.51.1.sif kallisto index -i hsGRCh38_ensv103_kallisto ${ENSV103_DIR:?unset}/Homo_sapiens.GRCh38.103.rna.fa.gz


```

## Transcript abundance quantification

The outputs with `kallisto quant` results for every sample are located on **figshare** dataset [Expression data and differential expression analysis of CDS2 targeted SW837 and OMM2.5 xenografted cells](https://figshare.com/articles/dataset/Expression_data_and_differential_expression_analysis_of_i_CDS2_i_targeted_SW837_and_OMM2_5_xenografted_cells/28016804) inside the `kallisto` folder which is equivalent to the results generated into `analysis/kallisto` by the steps shown below.

### Get FASTQ files from the Mapped BAM files

Human BAM files generated with `STAR` for SW837 samples and cleaned Human BAM files from OMM2.5 xengrafts samples, were transformed to FASTQ for the transcript quantification. This was done using the script `bamtofastq_from_manifest.R` 

**INPUT**: 
- `metadata/manifests/7687_cram_manifest_INFO_from_iRODS_all_wbam_counts_qc_wNOD_bams_finalbams.txt`: Manifest file with the BAM files used to get FASTQ files from for the quantification

```bash
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq
STUDY=7687
METADATADIR=${PROJECTDIR}/metadata
cd ${PROJECTDIR}

# Load environment with requring 
source ${PROJECTDIR}/scripts/expression/source_me.sh

# Path to the human reference genome used
HUMSQE_REF="/lustre/scratch125/core/sciops_repository/references/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa"

#This script takes the manifest 
Rscript ${PROJECTDIR}/scripts/expression/bamtofastq_from_manifest.R --manifest ${METADATADIR}/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_all_wbam_counts_qc_wNOD_bams_finalbams.txt --projectdir ${PROJECTDIR} --studyID ${STUDY:?unset} --mem 16000 --ref ${HUMSQE_REF:?unset}

```
**OUTPUTS**:
 - [`scripts/7687_bamtofastq_final_jobs.sh`](../scripts/7687_bamtofastq_final_jobs.sh) shell script with the commands to run the conversion of the BAM files to FASTQ files

#### Run the conversion of the BAM files to FASTQ files

```bash 
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq
STUDY=7687

cd ${PROJECTDIR}

# Load environment with requring 
source ${PROJECTDIR}/scripts/expression/source_me.sh

bash ${PROJECTDIR}/scripts/${STUDY}_bamtofastq_final_jobs.sh

#Place symbolic links to the main fastqs
cd ${PROJECTDIR}/fastqs/finalfastqs
find ../* -name "SW*" -exec cp -s {} ./ \;

```
**OUTPUTS**:
- FASTQ files for the quantification in the `fastqs/finalfastqs` directory

#### Generate the commands to run transcript abundance quantification with kallisto quant

The quantification was performed using the `kallisto quant` function. The script `run_kallisto_quant.sh` was used to run the quantification for all the samples.

**INPUT**:
- `metadata/manifests/7687_cram_manifest_INFO_from_iRODS_all_wbam_counts_qc_wNOD_bams_finalbams_fqs.txt`: manifest with the paths to the FASTQ files to use to quantify and sample names to be used for the output directories

```bash
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq
STUDY=7687
METADATADIR=${PROJECTDIR}/metadata
cd ${PROJECTDIR}

# Load environment with requring 
source ${PROJECTDIR}/scripts/expression/source_me.sh

## Kallisto singularity image location
KALLISTO_DIR=${PROJECTDIR}/scripts/singularity_images/kallisto
## Kalisto index location
KALLISTO_INDEX=${PROJECTDIR}/reference/kallisto/GRCh38_ensv103/hsGRCh38_ensv103_kallisto

# Load the modules to call singularity to generate the image 
module purge
module load singularity/3.11.4

#This script takes the cram manifest 
Rscript ${PROJECTDIR}/scripts/expression/make_kallisto_jobs_from_manifest.R --manifest ${METADATADIR:?unset}/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_all_wbam_counts_qc_wNOD_bams_finalbams_fqs.txt --projectdir ${PROJECTDIR:?unset} --studyID ${STUDY:?unset} --index ${KALLISTO_INDEX:?unset} --kallisto ${KALLISTO_DIR:?unset}/kallisto_0.51.1.sif --mem 16000

```
**OUTPUTS**:
- `scripts/7687_kallisto_final_jobs.sh` shell script with the `lsf bsub` commands to run the quantification of the samples

#### Run transcript abundance quantification with Kallisto quant

Singularity should be available in the running path to run the quantification 

```bash
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq
STUDY=7687
METADATADIR=${PROJECTDIR}/metadata
cd ${PROJECTDIR}

# Load environment with requring 
source ${PROJECTDIR}/scripts/expression/source_me.sh

## Kallisto singularity image location
KALLISTO_DIR=${PROJECTDIR}/scripts/singularity_images/kallisto
## Kalisto index location
KALLISTO_INDEX=${PROJECTDIR}/reference/kallisto/GRCh38_ensv103/hsGRCh38_ensv103_kallisto

# Load the modules to call singularity to generate the image 
module purge
module load singularity/3.11.4
# Runt he script
bash ${PROJECTDIR}/scripts/${STUDY}_kallisto_final_jobs.sh

```
**OUTPUTS**:
- The quantification results are located in the `analysis/kallisto`, one folder per sample is generated containing the following files : 
  - `abundance.h5` and `abundance.tsv` files for each sample
  - `run_info.json` file with the information of the run

- After all the jobs completed `abundance.tsv` files were compressed, the commands used were the follwing:

```bash
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq
STUDY=7687
cd ${PROJECTDIR}

# Load environment with requring 
source ${PROJECTDIR}/scripts/expression/source_me.sh

## KAllisto INDEX
KALLISTO_RES_DIR=${PROJECTDIR}/analysis/kallisto

# gzip all the files 
find  ${KALLISTO_RES_DIR:?unset}/*/* -name "abundance.tsv" -exec gzip {} \;
```
[!NOTE]  
:warning: **IMPORTANT NOTE** :warning:
Due to the file sizes of the quantification results, the outputs with `kallisto quant` results for every sample are located on **figshare** dataset [Expression data and differential expression analysis of CDS2 targeted SW837 and OMM2.5 xenografted cells](https://figshare.com/articles/dataset/Expression_data_and_differential_expression_analysis_of_i_CDS2_i_targeted_SW837_and_OMM2_5_xenografted_cells/28016804) inside the `kallisto` folder which is equivalent to the results generated into `analysis/kallisto` by the steps shown above.



