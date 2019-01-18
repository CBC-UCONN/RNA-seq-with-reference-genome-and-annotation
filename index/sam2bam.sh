#!/bin/bash
#SBATCH --job-name=sam2bam
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`


#################################################################
# SAM to BAM format
#################################################################
module load samtools/1.7

samtools view -@ 8 -uhS trimmed_LB2A_SRR1964642.sam | samtools sort -@ 4 -o sort_trim_LB2A_SRR1964642.bam
samtools view -@ 8 -uhS trimmed_LB2A_SRR1964643.sam | samtools sort -@ 4 -o sort_trim_LB2A_SRR1964643.bam
samtools view -@ 8 -uhS trimmed_LC2A_SRR1964644.sam | samtools sort -@ 4 -o sort_trim_LC2A_SRR1964644.bam
samtools view -@ 8 -uhS trimmed_LC2A_SRR1964645.sam | samtools sort -@ 4 -o sort_trim_LC2A_SRR1964645.bam


