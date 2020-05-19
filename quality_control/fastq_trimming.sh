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
# Trimming of Reads using Trimmomatic 
#################################################################
module load Trimmomatic/0.39


java -jar $Trimmomatic SE
	-threads 12 \
	../raw_data/LB2A_SRR1964642.fastq.gz \
	LB2A_SRR1964642_trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
	SLIDINGWINDOW:4:20 \
	MINLEN:45

java -jar $Trimmomatic SE
	-threads 12 \
	../raw_data/LB2A_SRR1964643.fastq.gz \
	LB2A_SRR1964643_trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
	SLIDINGWINDOW:4:20 \
	MINLEN:45

java -jar $Trimmomatic SE
	-threads 12 \
	../raw_data/LC2A_SRR1964644.fastq.gz \
	LC2A_SRR1964644_trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
	SLIDINGWINDOW:4:20 \
	MINLEN:45

java -jar $Trimmomatic SE
	-threads 12 \
	../raw_data/LC2A_SRR1964644.fastq.gz \
	LC2A_SRR1964644_trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
	SLIDINGWINDOW:4:20 \
	MINLEN:45


