#!/usr/bin/env bash

# Install miniconda

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
./miniconda3/bin/conda init bash

### ATTENTION!!! Here you should restart bash session

# Activate conda

conda create --yes --name papolarityenv python=3.7;
conda activate papolarityenv;

# Install necessary software conda

conda install --yes -c conda-forge parallel;
conda install --yes -c bioconda star bedtools samtools csvtk sra-tools;
pip install cutadapt pasio papolarity;
