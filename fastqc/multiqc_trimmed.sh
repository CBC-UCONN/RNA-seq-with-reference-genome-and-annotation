#!/bin/bash
#SBATCH --job-name=multiqc_trimmed
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
# MULTIQC of trimmed reads 
#################################################################
module load MultiQC/1.1

multiqc --outdir trimmed_multiqc ./after/