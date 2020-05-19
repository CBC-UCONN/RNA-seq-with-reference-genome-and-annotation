#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# FASTQC of raw reads 
#################################################################
dir="before"
if [ ! -d "$dir" ]; then
        mkdir -p $dir
fi

module load fastqc/0.11.5
fastqc --outdir ./"$dir"/ ../raw_data/LB2A_SRR1964642.fastq.gz
fastqc --outdir ./"$dir"/ ../raw_data/LB2A_SRR1964643.fastq.gz
fastqc --outdir ./"$dir"/ ../raw_data/LC2A_SRR1964644.fastq.gz
fastqc --outdir ./"$dir"/ ../raw_data/LC2A_SRR1964645.fastq.gz

