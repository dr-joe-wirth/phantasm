
# PHANTASM: <ins>PH</ins>ylogenomic <ins>AN</ins>alyses for the <ins>TA</ins>xonomy and <ins>S</ins>ystematics of <ins>M</ins>icrobes

## A docker image for this software can be found at https://hub.docker.com/r/jwirth/phantasm

## If you use our software, please cite our paper in Nucleic Acids Research
**Automating microbial taxonomy workflows with PHANTASM: PHylogenomic ANalyses for the TAxonomy and Systematics of Microbes**

Joseph S. Wirth & Eliot C. Bush, 2023

doi: [10.1093/nar/gkad196](https://doi.org/10.1093/nar/gkad196)

##### Note: PHANTASM requires an email address for communicating with NCBI as described in python's [Bio.Entrez](https://biopython.org/docs/latest/api/Bio.Entrez.html) package. The email address is not stored or used for any other purposes.

[![DOI](https://zenodo.org/badge/457891194.svg)](https://zenodo.org/badge/latestdoi/457891194)

## Table of Contents
1. [Installing PHANTASM](#1-installing-phantasm)

    1.1. [Using the Docker image](#11-using-the-docker-image)
    
    1.2. [Native Installation](#12-native-installation)
    
2. [Running PHANTASM](#2-running-phantasm)

3. [Dependencies](#3-dependencies)

    3.1. [External software](#31-external-software)
    
    3.2. [Python modules](#32-python-modules-python39-or-above)
    
    3.3 [R packages](#33-r-packages-r-version-411-or-above)

## 1. Installing PHANTASM
### 1.1 Using the Docker image
The easiest way to get started is is to use the [Docker image](https://hub.docker.com/r/jwirth/phantasm). There are instructions for using phantasm as a Docker container on [dockerhub](https://hub.docker.com/r/jwirth/phantasm). The Docker image contains all of the dependencies pre-installed and only requires that Docker Desktop is installed (or `docker` if using a linux server). Alternatively, Singularity may be used. Please see the [detailed tutorial for using a PHANTASM in a container](https://github.com/dr-joe-wirth/phantasm/blob/master/docker_build_files/README.md) for more information.


### 1.2 Native installation
PHANTASM requires `python3` (v. 3.9+). This can be installed via the command line with the following command:

    sudo apt-get -y install python3.9 python3-pip python3-setuptools python3-dev
          
`mummer` is needed for calculating ANI using [pyani](https://github.com/widdowquinn/pyani). It can be installed with the following command:

    sudo apt-get -y install mummer

`curl` is necessary for certain R libraries to function properly. It can be installed with the following command:

    sudo apt-get -y install curl libcurl4-openssl-dev

`blast+` is required for calculating core genes and average amino acid identity. It can be installed with the following command:

    sudo apt-get -y install ncbi-blast+

`R` needs to be installed for plotting heatmaps of the overall genome-related indices. It can be installed with the following command:

    sudo apt-get -y install r-base

The following R libraries need to be installed: `ape`, `gplots`, `dendextend`, and `BiocManager`. It also needs the `DECIPHER` bioconductor package. These can be installed using the following commands. _**Warning**: These packages can take several minutes to install._

    R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "install.packages('gplots',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "install.packages('dendextend',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
    R -e "BiocManager::install('DECIPHER')"

The following python packages need to be installed: `numpy`, `scipy`, `Bio`, `textdistance`, `parasail`, `rpy2`, `pyani`. They can be installed with `pip3` using the following commands:

    pip3 install scipy
    pip3 install Bio
    pip3 install textdistance
    pip3 install numpy
    pip3 install parasail
    pip3 install rpy2
    pip3 install pyani

PHANTASM requires [FastTree](https://www.microbesonline.org/fasttree/) (or FastTreeMP), [muscle v5.1]([MUSCLE v5.1](https://github.com/rcedgar/muscle/releases/tag/5.1.0)), and [iqtree v1.6.12](http://www.iqtree.org/#download) to be installed. It also needs the source code for [xenoGI v3.1.2](https://github.com/dr-joe-wirth/xenoGI/releases/tag/v3.1.2).

Finally, clone the PHANTASM repository or download a release.

In order for PHANTASM to function properly, you will need to modify the file `param.py` found in the `phantasm` directory. Specifically, you need to modify the following fields as shown in the table below:

     _______________________________________________________________________
    | variable name |                    variable value                     |
    |---------------|-------------------------------------------------------|
    | BLASTPLUS_DIR | absolute path to the blast+ executable directory (bin)|
    | MUSCLE_EXE    | absolute path the the muscle executable file          |
    | FASTTREE_EXE  | absolute path to the fasttree executable file         |
    | IQTREE_EXE    | absolute path to the IQTree executable file           |
    | PHANTASM_DIR  | absolute path to the phantasm directory               |
    | XENOGI_DIR    | absolute path to the xenoGI directory                 |
    |_______________________________________________________________________|

There are also optional values in `param.py` that can be modified in order to customize how PHANTASM runs. See the [docker tutorial](https://github.com/dr-joe-wirth/phantasm/blob/master/docker_build_files/README.md#22-modifying-phantasms-settings-optional) for more information.

To check that dependencies are properly installed, run the following command:

    python3 <path to phantasm>/phantasm.py check


## 2. Running PHANTASM
A [detailed tutorial](https://github.com/dr-joe-wirth/phantasm/blob/master/docker_build_files/README.md#3-running-phantasm) exists for the Docker image. The only major difference between running phantasm natively and running phantasm in a Docker container is that the command `phantasm` should be replaced with `python3 <path to phantasm>/phantasm.py`. Alternatively, you can add the following line to your bash profile (.profile, .bashrc, etc.) which will allow your terminal to call phantasm with the command `phantasm`. Keep in mind that you will need to reload your terminal for the changes to take effect.

    alias phantasm='python3 <absolute path to phantasm>/phantasm.py'

A help message can be obtained with the command `python3 <path to phantasm>/phantasm.py help`. Help messages for a given task can be obtained with the command `python3 <path to phantasm>/phantasm.py <TASK> --help`. The version of PHANTASM can be obtained with the command `python3 <path to phantasm>/phantasm.py version`. For help on analyzing the results, please see the [detailed tutorial](https://github.com/dr-joe-wirth/phantasm/blob/master/docker_build_files/README.md#4-analyzing-the-results).

## 3. Dependencies
See [section 1.2](#12-native-installation) for installation details.

### 3.1 External software
  * [FastTree](http://www.microbesonline.org/fasttree/)
  * [MUSCLE v5.1](https://github.com/rcedgar/muscle/releases/tag/5.1.0)
  * [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  * [IQTree](http://www.iqtree.org/#download)
  * [mummer](https://github.com/mummer4/mummer)

### 3.2 Python modules (python3.9 or above)
  * [xenoGI v3.1.2](https://github.com/dr-joe-wirth/xenoGI/releases/tag/v3.1.2)
  * scipy
  * Bio
  * textdistance
  * numpy
  * parasail
  * rpy2
  * [pyani](https://github.com/widdowquinn/pyani)
  * semver

### 3.3 R packages (R version 4.1.1 or above)
  * ape
  * gplots
  * dendextend
  * DECIPHER
