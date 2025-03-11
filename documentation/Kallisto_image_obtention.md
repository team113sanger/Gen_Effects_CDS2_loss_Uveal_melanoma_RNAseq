# Expression assessment using Kallisto

## Introduction

This document describes the steps to obtain a singularity image with the software Kallisto. This image will be used to run the expression assessment of the RNAseq data from the project **Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq**.

### Image generation

To do this we generated a singularity image using singularity [`v3.11.4`](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#installation)from the docker container which has [kallisto v0.51.1](https://github.com/pachterlab/kallisto/tree/master) from the biocontainers repository.

#### Instructions to generate the image

To install singularity `v3.11.4`, you can follow the instructions in the [official documentation](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#installation).

Then to generate the image, you can run the following commands:

```bash
PROJECTDIR=/lustre/7687_3364_Gen_Effects_CDS2_loss_Uveal_melanoma_RNAseq

cd ${PROJECTDIR}

#the docker  hub 
KALLISTO_DIR=${PROJECTDIR}/scripts/singularity_images/kallisto
mkdir -p ${KALLISTO_DIR}
cd  ${KALLISTO_DIR}

# Load the modules to call singularity to generate the image 
#module purge
#module load singularity/3.11.4

singularity pull --name kallisto_0.51.1.sif docker://quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0

# Test running command with the image
singularity exec -e -B `pwd` -B / -B ${PROJECTDIR:?unset}  \
${KALLISTO_DIR:?unset}/kallisto_0.51.1.sif kallisto --h

```