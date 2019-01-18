#!/bin/bash
#SBATCH --job-name=fastq_trimming
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Trimming of Reads using Sickle 
#################################################################
module load sickle/1.33 

sickle se -f ../raw_data/LB2A_SRR1964642.fastq -t sanger -o trimmed_LB2A_SRR1964642.fastq -q 30 -l 50
sickle se -f ../raw_data/LB2A_SRR1964643.fastq -t sanger -o trimmed_LB2A_SRR1964643.fastq -q 30 -l 50
sickle se -f ../raw_data/LC2A_SRR1964644.fastq -t sanger -o trimmed_LC2A_SRR1964644.fastq -q 30 -l 50 
sickle se -f ../raw_data/LC2A_SRR1964645.fastq -t sanger -o trimmed_LC2A_SRR1964645.fastq -q 30 -l 50
