#!/bin/bash
#SBATCH --job-name=fastqc_trimmed
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
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
dir="after"
if [ ! -d "$dir" ]; then
        mkdir -p $dir
fi

module load fastqc/0.11.5
fastqc --outdir ./"$dir"/ ../quality_control/trimmed_LB2A_SRR1964642.trim.fastq.gz -t 8
fastqc --outdir ./"$dir"/ ../quality_control/trimmed_LB2A_SRR1964643.trim.fastq.gz -t 8
fastqc --outdir ./"$dir"/ ../quality_control/trimmed_LC2A_SRR1964644.trim.fastq.gz -t 8
fastqc --outdir ./"$dir"/ ../quality_control/trimmed_LC2A_SRR1964645.trim.fastq.gz -t 8

