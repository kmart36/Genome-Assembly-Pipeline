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

Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." Cell Systems 3(1), 2016.

Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021) Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nat Methods, 18:170-175. https://doi.org/10.1038/s41592-020-01056-5

Cheng, H., Jarvis, E.D., Fedrigo, O., Koepfli, K.P., Urban, L., Gemmell, N.J., Li, H. (2022) Haplotype-resolved assembly of diploid genomes without parental data. Nature Biotechnology, 40:1332–1335. https://doi.org/10.1038/s41587-022-01261-x

Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]

Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654

Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323

Guillaume Marcais and Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics (2011) 27(6): 764-770 (first published online January 7, 2011) doi:10.1093/bioinformatics/btr011