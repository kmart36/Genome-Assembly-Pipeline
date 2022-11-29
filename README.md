# Running the Genome Assembly Pipeline
## Installation
This pipeline is designed to be as simple to use as possible: input your data, and a completed, high quality assembly will come out. (MORE INTRO) This tutorial will walk through how to install and use every program listed. Links to the programs in the pipeline are below. 

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
