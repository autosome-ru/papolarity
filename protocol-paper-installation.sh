#!/usr/bin/env bash

# Install miniconda

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh  # At this stage some manual actions from the user are required.

# if you copy-paste commands don't copy the following commands into the shell before the previous one is completed
./miniconda3/bin/conda init bash


### ATTENTION!!! Here you should restart bash session

# Activate conda

conda create --yes --name papolarityenv python=3.7;
conda activate papolarityenv;

# Install necessary software conda

conda install --yes -c conda-forge parallel=20191122;
conda install --yes -c bioconda star=2.7.3a bedtools=2.29.2 samtools=1.9 csvtk=0.19.1 sra-tools=2.8.0;
pip install cutadapt==2.7 pasio==1.1.2 papolarity==1.0.0;
