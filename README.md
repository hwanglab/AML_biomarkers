[![Open in Visual Studio Code](https://open.vscode.dev/badges/open-in-vscode.svg)](https://open.vscode.dev/hwanglab/AML_biomarkers)

# Biomarker Analysis for Single-Cell Data

## Installation/Setup

Note: You also need to have access to rsinglecell--a custom package. 

### Using Docker Image
The simplest way to run this with Docker is using Visual Studio Code--but most of these steps work without VSCode.

Prior to building the container, you must do 2 things.
Firstly, create a file called `.env` in the root directory and add your GitHub email (VAR:`EMAIL`), PAT (VAR:`GIT_TOKEN`), your name (VAR:`NAME`), and your GitHub username (VAR:`UNAME`).
Your email should also be the same as used with CIBERSORTx.
You should also add your CIBERSORTx token (VAR:`TOKEN`).
This file will be automatically recognized at build time.
Secondly, edit the `docker-compose.yml` file with the following line depending on your host system (do not use realitive paths):
```
# Windows
# %LOCALAPPDATA% is C:\Users\[Your Username]\AppData\Local
%LOCALAPPDATA%/renv:/root/.local/share/renv:cached

# MacOS (notice the space in Application Support)
/Users/[Your Username]/Library/Application Support/renv:/root/.local/share/renv:cached

# Linux
/usr/[Your Username]/.local/share/renv:/root/.local/share/renv:cached

# If you have modified $RENV_CACHE_ROOT on any system
$RENV_CACHE_ROOT:/root/.local/share/renv:cached
```
This will **not** effect your local R packages--even if they are of a different R version, or system.
If you do not use [renv](https://rstudio.github.io/renv), please create this folder.
You can do so by using renv, or via some folder wizardry.
When starting VSCode with the Remote Containers extension, build the image. 
This will take some time as all the R packages will be installed.
From there, you should have an interactive development enviorment and can run the pipeline. 
Due to caching by renv, successive builds should be faster.
And, if you build other containers like this, and use the same packages, they too will build faster.

#### Running without Visual Studio Code
Using VSCode makes running simpler, as directory binding is somewhat automatic. 
If you decide against using VSCode, you first need to uncomment the following line in the Dockerfile: `RUN Rscript -e "source('renv/activate.R');renv::restore(prompt = FALSE)"`.
Then, you have to take care to bind all the directories from your system into the container.
To mount a volume to use for the renv cache use the [approparte cache location](https://rstudio.github.io/renv/reference/paths.html) and bind with the following structure: `"$RENV_ROOT_HOST":/root/.local/share/renv`. 
Please note, the quotes are required if spaces are present in the path (Thanks Apple!). Then just use `docker compose up` to start the container.

### Getting Setup without Docker
Renv is used for this project. 
To install all the R packages just use `renv::restore()`. 
(This is what happens during the docker container build.)
The renv autoloader should automatically boostrap renv and `renv::restore()` will take care of installing all the correct package versions, although you will have to install system dependencies yourself.

### Getting set up on OSC
Singlularity was not trivial to use, therfore I directly loaded modules. The following should be added to your `.bashrc` file or prepend to any submission script:
```
# Set up R enviorment
export PATH=$PATH:$HOME/.local/bin #if you like to use radian
module load hdf5-serial/1.12.0
module load magick/version
module load fftw3/3.3.8
module load R/4.1.0
module load gdal/3.2.2
module load proj/7.2.1
module load geos/3.8.2
module load python/3.7-2019.10

# CIBERSORT options
export TOKEN={Token from website}
export EMAIL={email address used on website}
```

Then follow the directions under "Getting Setup without Docker"

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

**CIBERSORTx.py:**
    Runs the CIBERSORTx algorithm on the data.

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