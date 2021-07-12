[![Open in Visual Studio Code](https://open.vscode.dev/badges/open-in-vscode.svg)](https://open.vscode.dev/hwanglab/AML_biomarkers)

# Biomarker Analysis for Single-Cell Data

## Installation/Setup

### Using Docker Image
The simplest way to run this with Docker is using Visual Studio Code.
When starting VSCode with the Remote Containers extension, build the image. 
This will take some time as all the R packages will be installed.
From there, you should have an interactive development enviorment and can run the pipeline. 
Due to caching by renv, successive builds should be faster.

#### Running without Visual Studio Code
Using VSCode makes running simpler, as directory binding is somewhat automatic. 
If you decide against using VSCode, you first need to uncomment the following line in the Dockerfile: `RUN Rscript -e "source('renv/activate.R');renv::restore(prompt = FALSE)"`.
Then, you have to take care to bind all the directories from your system into the container.

### Getting Setup without Docker
Renv is used for this project. 
To install all the R packages just use `renv::restore()`. 
(This is what happens during the docker container build.)
The renv autoloader should automatically boostrap renv and `renv::restore()` will take care of installing all the correct package versions, although you will have to install system dependencies yourself.

Note: You also need to have access to rsinglecell--a custom package. 

## Running the Pipeline

### General Overview
The pipeline is stored in the cli directory. 
The functions each take an id which will automatically handle inputs and outputs. 
An invocation file is created if not present and will contain the arguments used at runtime. 
CLI functions cache some of their outputs to make repeated runs quicker. 
All outputs are saved, however, when cached, the expressions are only rerun if the data changes.
To make the deconvolution references, upload the cibersort_ref_input.txt file to CIBERSORTx and then download the inferred_ref file from CIBERSORTx.
This file should start with CIBERSORTx_Job and be in the output directory for the pipeline instance.

### Data Required
Some external datasets are required. TARGET AML, BeatAML, and TCGA LAML are all used during deconvolution. 
Many different clinical information tables are also required. 
These are all found in `/clinical_info`

### CLI Function Information

#### Function Descriptions
**dimension_reduction.R:**
    Subsets the data and runs dimension reductions (UMAP, PCA) and does batch correction using Harmony. 
    Then finds clusters that are different between groups.
    The results are cached.

**make_CIBERSORT_ref.R:**
    Prepares data to upload to CIBERSORTx for reference creation. 
    The results are cached.

**deconvolute.R:**
    Prepares references for deconvolution and runs several deconvolution methods.
    The results are cached.

**prepare_deconvoluted_samples.R:**
    Annotates deconvoluted samples with clinical information. Selects desired deconvolution method. 

**run_gsva.R:**
    Runs a GSVA on the clusters.
    The results are cached.

**run_cellphone.R:**
    Runs the CellPhoneDB algorithm on the clusters.
    The results are cached.

**survival_analysis.R:**
    Does a survival analysis using the chosen parameters.

#### Common CLI Arguments

`--id`: An arbitrary ID to give the run. 
Should be consistant across all functions. 
Will be a subfolder of `--dir`.

`--dir`: What directory to use as the run directory. 
Defualts to the current directory.

`--verbose`: How should messages be printed.

`--cores`: How many cores should be used for multicore processing. 
This uses the future framework.

`--help`: will print help for each CLI function.  

## Future Goals
- Simplify Clinical Annotation sheets
- Clear Downloading of External Datasets