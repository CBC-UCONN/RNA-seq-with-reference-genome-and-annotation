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
module load hisat2/2.2.0

hisat2 -p 8 --dta -x ../index/L_crocea -U ../quality_control/LB2A_SRR1964642_trim.fastq.gz | \
	samtools view -S -h -u - | \
	samtools sort -T SRR1964642 - >LB2A_SRR1964642.bam

hisat2 -p 8 --dta -x ../index/L_crocea -U ../quality_control/LB2A_SRR1964643_trim.fastq.gz | \
	samtools view -S -h -u - | \
	samtools sort -T SRR1964643 - >LB2A_SRR1964643.bam

hisat2 -p 8 --dta -x ../index/L_crocea -U ../quality_control/LC2A_SRR1964644_trim.fastq.gz | \
	samtools view -S -h -u - | \
	samtools sort -T SRR1964644 - >LC2A_SRR1964644.bam

hisat2 -p 8 --dta -x ../index/L_crocea -U ../quality_control/LC2A_SRR1964645_trim.fastq.gz | \
	samtools view -S -h -u - | \
	samtools sort -T SRR1964645 - >LC2A_SRR1964645.bam

