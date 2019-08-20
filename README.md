[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3361700.svg)](https://doi.org/10.5281/zenodo.3361700)
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.0-brightgreen.svg)](https://snakemake.bitbucket.io)

# Data analysis for the initial Varlociraptor paper

This Snakemake workflow automatically generates all results and figures from the paper.

## Requirements

Any 64-bit Linux installation with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6).
Note that the restriction of this workflow to Linux is purely a design decision (to save space and ensure reproducibility) and not related to Conda/Bioconda. Bioconda packages are available for both Linux and MacOS in general.

In this special case there are unfortunately unavoidable additional requirements, due to licensing restrictions and data size.

1. The required real data has to be obtained from EGA: https://ega-archive.org/datasets/EGAD00001002142
2. The [Lancet](https://github.com/nygenome/lancet) variant caller has to be manually installed (and available in your PATH). It cannot be automatically deployed by Snakemake due to licensing restrictions.
3. The required simulated data has to be obtained from Zenodo. Due to technical issues, we have not yet managed to upload it. This is being worked on with Zenodo support.

## Usage

This workflow can be used to recreate all results found in the paper.

### Step 1: Setup system

#### Variant a: Installing Miniconda on your system

If you are on a Linux system with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6), you can simply install Miniconda3 with

    curl -o /tmp/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash /tmp/miniconda.sh

Make sure to answer `yes` to the question whether your PATH variable shall be modified.
Afterwards, open a new shell/terminal.

#### Variant b: Use a Docker container

Otherwise, e.g., on MacOS or if you don't want to modify your system setup, install [Docker](https://www.docker.com/), run

    docker run -it continuumio/miniconda3 /bin/bash
  
and execute all the following steps within that container.

#### Variant c: Use an existing Miniconda installation

If you want to use an existing Miniconda installation, please be aware that this is only possible if it uses Python 3 by default. You can check this via
  
    python --version

Further, ensure it is up to date with

    conda update --all

### Step 2: Setup Bioconda channel

Setup Bioconda with

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

### Step 3: Install Snakemake

Install Snakemake >=5.4.0 with

    conda install snakemake

If you already have an older version of Snakemake, please make sure it is updated to >=5.4.0.

### Step 4: Download the workflow

First, create a working directory:

    mkdir varlociraptor-workflow
    cd varlociraptor-workflow

Then, download the workflow archive from https://doi.org/10.5281/zenodo.3361700 and unpack it with

    tar -xf workflow.tar.gz

### Step 5: Run the workflow

Execute the analysis workflow with Snakemake

    snakemake --use-conda

Please wait a few minutes for the analysis to finish.
Results can be found in the folder `figs/`.
If you have been running the workflow in the docker container (see above), 
you can obtain the results with

    docker cp <container-id>:/bioconda-workflow/figs .

whith `<container-id>` being the ID of the container.


## Known errors

* If you see an error like
  ```
  ImportError: No module named 'appdirs'
  ```
  when starting Snakemake, you are likely suffering from a bug in an older conda version. Make sure to update your conda installation with 

      conda update --all

  and then reinstall the `appdirs` and `snakemake` package with

      conda install -f appdirs snakemake
* If you see an error like
  ```
  ImportError: Missing required dependencies ['numpy']
  ```
  you are likely suffering from a bug in an older conda version. Make sure to update your conda installation with
  
      conda update --all
  
  and then reinstall the `snakemake` package with

      conda install -f snakemake
