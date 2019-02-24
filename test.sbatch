#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH -p RM
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH --output=test.log
#SBATCH --chdir=/home/siweixu/project

source ~/project/src/load_module.sh
source main.sh -l ~project/data/test_1.fastq -r ~/project/data/test_2.fastq -g ~/project/data/genome/fasta/hg38.fa -c ~/project/data/genome.gtf.gz
