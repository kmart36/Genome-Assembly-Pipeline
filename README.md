# Running the Genome Assembly Pipeline
## Installation
This pipeline is designed to be as simple to use as possible: input your data, and a completed, high quality assembly will come out. Currently, this pipeline is only useable on Bucknell University's BisonNet. It also utilizes the SLURM system which allows programs to run on a system with much more computing power than a local computer. This tutorial will walk through how to install and use every program listed. Links to the programs in the pipeline are below. 

#### Programs Used:
* Nextflow: https://www.nextflow.io/docs/latest/getstarted.html 
* Miniconda: https://docs.conda.io/en/latest/miniconda.html 
* fasterq-dump: https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump 
* fastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
* jellyfish: https://github.com/gmarcais/Jellyfish 
* hifiasm: https://github.com/chhylp123/hifiasm 
* BUSCO: https://busco.ezlab.org/ 
* bwa: https://bio-bwa.sourceforge.net/bwa.shtml 
* SALSA2: https://github.com/marbl/SALSA 
* Juicer: https://github.com/aidenlab/juicer/wiki/Download 

Many scripts will be provided for you. This guide will tell you step-by-step what you will need to edit, if necessary, and how to run each program. 

#### Getting started: 

Nextflow is key to running the entire pipeline. Nextflow is a bioinformatics workflow manager that allows for parallelization of programs as well as the simplicity of streaming together processes. 

Make sure you are in your “outermost” directory. For those on the BisonNet system, this will be the directory that looks like your Bucknell login. If you would like to keep this process separate from your other work, you can run these commands to create a separate directory. 
```
mkdir genome-assembly
cd genome-assembly
```
You can change `genome-assembly` to whatever you would like the folder to be called. The `cd` command ensures that you are working inside your new folder. 

For downloading Nextflow, run these two commands in your Linux terminal. 
```
wget -q0- https://get.nextflow.io | bash
chmod +x nextflow
```
When you are finished, your current directory should have a file called `nextflow`. This is the main executable file and will allow you to run your Nextflow scripts. 

Next, to download the script, click here (eventually will publish to github).
(Insert instructions)

For most of our processes, we use miniconda which allows for easy installs and easy usage. To install from the command line, run these commands in your Linux terminal. 
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
#### The following programs will need to be installed through miniconda:
* hifiasm
* bwa
* assembly-stats

#### The following programs will need to be installed separately: 
* SALSA2
* Juicer

#### The following programs are already installed on the BisonNet system in the form of modules:
* fasterq-dump
* fastQC
* Jellyfish
* BUSCO
