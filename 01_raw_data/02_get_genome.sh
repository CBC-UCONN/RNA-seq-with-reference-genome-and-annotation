#!/bin/bash
#SBATCH --job-name=get_genome
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

#################################################################
# Download genome and annotation from ENSEMBL
#################################################################

# load software
module load samtools/1.12

# output directory
GENOMEDIR=../genome
mkdir -p $GENOMEDIR

# we're using Fundulus heteroclitus from ensembl v105
    # we'll download the genome, GTF annotation and transcript fasta
    # https://useast.ensembl.org/Fundulus_heteroclitus/Info/Index


# download the genome
wget ftp://ftp.ensembl.org/pub/release-105/fasta/fundulus_heteroclitus/dna/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.dna.toplevel.fa.gz
# decompress it
gunzip Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.dna.toplevel.fa.gz

# download the GTF annotation
wget ftp://ftp.ensembl.org/pub/release-105/gtf/fundulus_heteroclitus/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.105.gtf.gz
# decompress it
gunzip Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.105.gtf.gz

# download the transcript fasta
wget http://ftp.ensembl.org/pub/release-105/fasta/fundulus_heteroclitus/cds/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cds.all.fa.gz
# decompress it
gunzip Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cds.all.fa.gz

# generate simple samtools fai indexes 
samtools faidx Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.dna.toplevel.fa
samtools faidx Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cds.all.fa

# move everything to the genome directory
mv Fundulus* $GENOMEDIR