#!/bin/bash
#SBATCH --job-name=htseq_count
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

echo `hostname`
date

#################################################################
# Generate Counts 
#################################################################
module load htseq/0.13.5
module load parallel/20180122

INDIR=../04_align/alignments
OUTDIR=counts
mkdir -p $OUTDIR

# accession list
ACCLIST=../01_raw_data/accessionlist.txt

# gtf formatted annotation file
GTF=../genome/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.105.gtf

# run htseq-count on each sample, up to 5 in parallel
cat $ACCLIST | \
parallel -j 5 \
    "htseq-count \
        -s no \
        -r pos \
        -f bam $INDIR/{}.bam \
        $GTF \
        > $OUTDIR/{}.counts"

