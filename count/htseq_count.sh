#!/bin/bash
#SBATCH --job-name=htseq_count
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
# Download gff
#################################################################
wget ftp://ftp.ensembl.org/pub/release-100/gtf/larimichthys_crocea/Larimichthys_crocea.L_crocea_2.0.100.gtf.gz

gunzip *.gz


#################################################################
# Generate Counts 
#################################################################
module load htseq/0.11.2

htseq-count -s no -r pos -f bam ../align/LB2A_SRR1964642.bam Larimichthys_crocea.L_crocea_2.0.100.gtf > LB2A_SRR1964642.counts
htseq-count -s no -r pos -f bam ../align/LB2A_SRR1964643.bam Larimichthys_crocea.L_crocea_2.0.100.gtf > LB2A_SRR1964643.counts
htseq-count -s no -r pos -f bam ../align/LC2A_SRR1964644.bam Larimichthys_crocea.L_crocea_2.0.100.gtf > LC2A_SRR1964644.counts
htseq-count -s no -r pos -f bam ../align/LC2A_SRR1964645.bam Larimichthys_crocea.L_crocea_2.0.100.gtf > LC2A_SRR1964645.counts
