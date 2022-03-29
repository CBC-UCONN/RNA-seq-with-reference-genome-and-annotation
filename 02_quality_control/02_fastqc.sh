#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastqc/0.11.7
module load parallel/20180122

# set input/output directory variables
INDIR=trimmed_sequences
REPORTDIR=fastqc_reports
mkdir -p $REPORTDIR

ACCLIST=../01_raw_data/accessionlist.txt
# run fastp in parallel, 4 samples at a time
cat $ACCLIST | parallel -j 4 \
    fastqc --outdir $REPORTDIR $INDIR/{}_trim_{1..2}.fastq.gz