#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`


#################################################################
# Aligning to Genome
#################################################################
module load hisat2/2.0.5

hisat2 -p 8 --dta -x ../index/L_crocea -q ../quality_control/trimmed_LB2A_SRR1964642.fastq -S trimmed_LB2A_SRR1964642.sam
hisat2 -p 8 --dta -x ../index/L_crocea -q ../quality_control/trimmed_LB2A_SRR1964643.fastq -S trimmed_LB2A_SRR1964643.sam
hisat2 -p 8 --dta -x ../index/L_crocea -q ../quality_control/trimmed_LC2A_SRR1964644.fastq -S trimmed_LC2A_SRR1964644.sam
hisat2 -p 8 --dta -x ../index/L_crocea -q ../quality_control/trimmed_LC2A_SRR1964645.fastq -S trimmed_LC2A_SRR1964645.sam
