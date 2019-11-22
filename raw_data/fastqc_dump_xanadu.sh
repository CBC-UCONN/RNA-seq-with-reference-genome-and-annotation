#!/bin/bash
#SBATCH --job-name=fastq_dump_xanadu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Download fastq files from SRA 
#################################################################
module load sratoolkit/2.8.2

fastq-dump SRR1964642
fastq-dump SRR1964643
fastq-dump SRR1964644
fastq-dump SRR1964645

#################################################################
# Rename the files 
#################################################################
mv SRR1964642.fastq LB2A_SRR1964642.fastq
mv SRR1964643.fastq LB2A_SRR1964643.fastq
mv SRR1964644.fastq LC2A_SRR1964644.fastq
mv SRR1964645.fastq LC2A_SRR1964645.fastq

