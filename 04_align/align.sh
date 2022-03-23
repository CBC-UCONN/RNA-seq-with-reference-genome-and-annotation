#!/bin/bash
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=[0-18]%5

echo `hostname`

#################################################################
# Align reads to genome
#################################################################
module load hisat2/2.2.1
module load samtools/1.12

INDIR=../02_quality_control/trimmed_sequences
OUTDIR=alignments
mkdir -p $OUTDIR

# this is an array job. 
	# one task will be spawned for each sample
	# for each task, we specify the sample as below
	# use the task ID to pull a single line, containing a single accession number from the accession list
	# then construct the file names in the call to hisat2 as below

ACCLIST=../01_raw_data/accessionlist.txt

NUM=$(expr ${SLURM_ARRAY_TASK_ID} + 1)

SAMPLE=$(sed -n ${NUM}p $ACCLIST)

INDEX=../genome/hisat2_index/Fhet

# run hisat2
hisat2 \
	-p 8 \
	-x $INDEX \
	-1 $INDIR/${SAMPLE}_trim_1.fastq.gz \
	-2 $INDIR/${SAMPLE}_trim_2.fastq.gz | \
samtools view -@ 8 -S -h -u - | \
samtools sort -@ 8 -T $SAMPLE - >$OUTDIR/$SAMPLE.bam

# index bam files
samtools index $OUTDIR/$SAMPLE.bam
