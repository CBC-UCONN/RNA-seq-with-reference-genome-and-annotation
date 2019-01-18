#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Download the Genome
#################################################################
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/972/845/GCF_000972845.1_L_crocea_1.0/GCF_000972845.1_L_crocea_1.0_genomic.fna.gz

gunzip *.gz

#################################################################
# Indexing the Genome
#################################################################
module load hisat2/2.0.5
hisat2-build -p 16 GCF_000972845.1_L_crocea_1.0_genomic.fna L_crocea
