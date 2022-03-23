#!/bin/bash
#SBATCH --job-name=stringtie_merge
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
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
# Merge stringtie assemblies
#################################################################

# load software
module load stringtie/2.1.5

# input/output variables
INDIR=transcripts
OUTDIR=merged_transcripts
mkdir -p $OUTDIR

# stringtie GTF list
ls $INDIR/*gtf >gtflist.txt

# ENSEMBL GTF annotation file
GTF=../genome/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.105.gtf

# run stringtie merge
stringtie --merge \
    -p 10 \
    -G $GTF \
    -o $OUTDIR/merged.gtf \
    gtflist.txt
