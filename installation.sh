#!/bin/bash 
#SBATCH -p medium # partition (queue) 
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code) 
#SBATCH -n 8 # number of cores 
#SBATCH --mem-per-cpu=65536 # memory per core 
#SBATCH --job-name="installation" # job name 
#SBATCH -o slurm.%N.%j.stdout.txt # STDOUT 
#SBATCH -e slurm.%N.%j.stderr.txt # STDERR 
#SBATCH --mail-user=USERID@bucknell.edu # address to email 
#SBATCH --mail-type=ALL # mail events (NONE, BEGIN, END, FAIL, ALL) 

#Installing Nextflow
mkdir genome-assembly
cd genome-assembly
wget -q0- https://get.nextflow.io | bash
chmod +x nextflow

#Installing miniconda
cd ..
mkdir miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3/miniconda.sh
bash miniconda3/miniconda.sh -b -u -p miniconda3
rm -rf miniconda3/miniconda.sh
miniconda3/bin/conda init bash
miniconda3/bin/conda init zsh

#Installing hifiasm
conda install -c bioconda hifiasm

#Installing bwa
conda install -c bioconda bwa

#Installing assembly-stats
conda install -c bioconda assembly-stats

#Setting up git for SALSA and Juicer installation
mkdir SALSA
cd SALSA
git init
git clone git@github.com:marbl/SALSA.git

#
cd ..
