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
# Generated Counts 
#################################################################
module load htseq/0.9.1

htseq-count -s no -r pos -t gene -i Dbxref -f bam ../align/sort_trim_LB2A_SRR1964642.bam GCF_000972845.1_L_crocea_1.0_genomic.gff > LB2A_SRR1964642.counts
htseq-count -s no -r pos -t gene -i Dbxref -f bam ../align/sort_trim_LB2A_SRR1964643.bam GCF_000972845.1_L_crocea_1.0_genomic.gff > LB2A_SRR1964643.counts
htseq-count -s no -r pos -t gene -i Dbxref -f bam ../align/sort_trim_LC2A_SRR1964644.bam GCF_000972845.1_L_crocea_1.0_genomic.gff > LC2A_SRR1964644.counts
htseq-count -s no -r pos -t gene -i Dbxref -f bam ../align/sort_trim_LC2A_SRR1964645.bam GCF_000972845.1_L_crocea_1.0_genomic.gff > LC2A_SRR1964645.counts
