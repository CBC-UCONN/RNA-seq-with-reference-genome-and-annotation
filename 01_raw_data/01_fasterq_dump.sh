#!/bin/bash
#SBATCH --job-name=fastqer_dump_xanadu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`
date

# load parallel module
module load parallel/20180122

#################################################################
# Download fastq files from SRA 
#################################################################
# The data are a subset (2 populations) from this study:
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156460

# you need to have sratoolkit installed
    # easiest to do this locally
    # https://github.com/ncbi/sra-tools

cat accessionlist.txt | parallel -j 2 fasterq-dump

#################################################################
# compress the files 
#################################################################

ls *fastq | parallel -j 12 gzip
