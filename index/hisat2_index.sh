#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Download the Genome
#################################################################

# download it
wget ftp://ftp.ensembl.org/pub/release-100/fasta/larimichthys_crocea/dna/Larimichthys_crocea.L_crocea_2.0.dna.toplevel.fa.gz

# decompress it
gunzip Larimichthys_crocea.L_crocea_2.0.dna.toplevel.fa.gz

#################################################################
# Indexing the Genome
#################################################################

module load hisat2/2.2.0

hisat2-build -p 16 Larimichthys_crocea.L_crocea_2.0.dna.toplevel.fa L_crocea
