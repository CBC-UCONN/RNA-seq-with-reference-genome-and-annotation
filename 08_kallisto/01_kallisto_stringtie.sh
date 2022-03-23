#!/bin/bash
#SBATCH --job-name=kallisto_stringtie
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Quantify transcript abundance with Kallisto
#################################################################

# load software
module load kallisto/0.46.1
module load parallel/20180122
module load gffread/0.12.7
module load samtools/1.12

# input/output variables

INDIR=../02_quality_control/trimmed_sequences
STRINGTIE=quant_stringtie
mkdir -p $STRINGTIE

# accession list
ACCLIST=../01_raw_data/accessionlist.txt

# genome fasta
GENOME=../genome/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.dna.toplevel.fa

# stringtie gtf file
GTF=../07_stringtie/merged_transcripts/merged.gtf

# quantify expression against stringtie transcripts -------------------------------------------

# extract stringtie transcript sequences from genome fasta file
gffread -w stringtie_transcripts.fa -g $GENOME $GTF

# index transcript fasta file
kallisto index -i stringtie_transcripts.idx stringtie_transcripts.fa

# quantify transcript abundance for each sample against stringtie transcripts, up to 5 at a time
cat $ACCLIST | \
parallel -j 5 \
    "kallisto quant \
        -i stringtie_transcripts.idx \
        -o $STRINGTIE/{} \
        -b 100 \
        <(zcat $INDIR/{}_trim_1.fastq.gz) \
        <(zcat $INDIR/{}_trim_2.fastq.gz)"

