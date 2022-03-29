#!/bin/bash
#SBATCH --job-name=fastp_trimming
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
date

#################################################################
# Trimming/QC of reads using fastp
#################################################################
module load fastp/0.23.2
module load parallel/20180122

# set input/output directory variables
INDIR=../01_raw_data
REPORTDIR=fastp_reports
mkdir -p $REPORTDIR
TRIMDIR=trimmed_sequences
mkdir -p $TRIMDIR

ACCLIST=$INDIR/accessionlist.txt

# run fastp in parallel, 4 samples at a time
cat $ACCLIST | parallel -j 4 \
fastp \
	--in1 $INDIR/{}_1.fastq.gz \
	--in2 $INDIR/{}_2.fastq.gz \
	--out1 $TRIMDIR/{}_trim_1.fastq.gz \
	--out2 $TRIMDIR/{}_trim_2.fastq.gz \
	--json $REPORTDIR/{}_fastp.json \
	--html $REPORTDIR/{}_fastp.html

