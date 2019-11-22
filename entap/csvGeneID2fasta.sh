#!/bin/bash
#SBATCH --job-name=csvGeneID2fasta
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

head -n 10 Croaker_DESeq2-results-with-normalized.csv > temp.csv

module load biopython/1.70

python csvGeneID2fasta.py temp.csv ProteinTable12197_229515.txt GCF_000972845.1_L_crocea_1.0_protein.faa > GeneID_proteinID.txt

