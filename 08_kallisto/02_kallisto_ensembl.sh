#!/bin/bash
#SBATCH --job-name=kallisto
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
ENSEMBL=quant_ensembl
mkdir -p $ENSEMBL

# accession list
ACCLIST=../01_raw_data/accessionlist.txt

# quantify expression against ensembl transcripts -------------------------------------------

# ensemble transcripts
CDS=../genome/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cds.all.fa

# index transcript fasta file 
kallisto index -i ensembl_transcripts.idx $CDS

# quantify transcript abundance for each sample against ensembl sequences, up to 5 at a time
cat $ACCLIST | \
parallel -j 5 \
    "kallisto quant \
        -i ensembl_transcripts.idx \
        -o $ENSEMBL/{} \
        -b 100 \
        <(zcat $INDIR/{}_trim_1.fastq.gz) \
        <(zcat $INDIR/{}_trim_2.fastq.gz)"

# extract transcript to gene mapping for downstream analysis

# ensemble GTF
GTF=../genome/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.105.gtf

# use regular expressions to pull out IDs
paste \
    <(awk '$3 ~ /transcript/' $GTF | grep -oP "(?<=gene_id \")[A-Z0-9]+") \
    <(awk '$3 ~ /transcript/' $GTF | grep -oP "(?<=transcript_id \")[A-Z0-9]+") \
>ensembl_gene2tx.txt