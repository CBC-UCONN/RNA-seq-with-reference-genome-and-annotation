#!/bin/bash
#SBATCH --job-name=stringtie_assemble
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Assemble transcripts with stringtie
#################################################################

# load software
module load stringtie/2.1.5
module load parallel/20180122

# input/output variables

INDIR=../04_align/alignments
OUTDIR=transcripts
mkdir -p 

# accession list
ACCLIST=../01_raw_data/accessionlist.txt

# ENSEMBL GTF annotation file
GTF=../genome/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.105.gtf

# run stringtie on all samples, up to 5 in parallel
cat $ACCLIST | \
parallel -j 5 \
    "stringtie \
        -o $OUTDIR/{}.gtf \
        -G $GTF \
        -p 2 \
        $INDIR/{}.bam"