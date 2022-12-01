#!/bin/bash 
#SBATCH -p medium # partition (queue) 
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code) 
#SBATCH -n 4 # number of cores 
#SBATCH --mem-per-cpu=32768 # memory per core 
#SBATCH --job-name="hifi-hic" # job name 
#SBATCH -o slurm.%N.%j.stdout.txt # STDOUT 
#SBATCH -e slurm.%N.%j.stderr.txt # STDERR 
#SBATCH --mail-user=kam071@bucknell.edu # address to email 
#SBATCH --mail-type=ALL # mail events (NONE, BEGIN, END, FAIL, ALL) 

cd SampleReads/Rhagonycha_fulva
hifiasm -l 0 -o ERR6606788_hic_nopurge.asm -t 32 --h1 ERR6054569_1.fastq --h2 ERR6054569_2.fastq ERR6606788.fastq  
