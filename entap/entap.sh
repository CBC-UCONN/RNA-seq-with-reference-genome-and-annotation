#!/bin/bash
#SBATCH --job-name=entap
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

module load EnTAP/0.9.0-beta
module load diamond/0.9.19

EnTAP --runP -i fasta_out.fasta -d /isg/shared/databases/Diamond/RefSeq/vertebrate_other.protein.faa.97.dmnd -d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd --ontology 0  --threads 16 


