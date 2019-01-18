#!/usr/bin/env python

##########################################################################################
## Coded by Neranjan V Perera (2018 Feb 25)						##
## Computational Biology Core 								##
## Institute for Systems Genomics 							##
## University of Connecticut 								##
## Copyright  										##
##											##
## The program will read the csv file with genes and give a fasta sequence 		##
## of the protein sequence								##
##											##
## usage: 										##
##	python csvGeneID2fasta.py <csv file>  <tabular file> <protein fasta>		##
##	csv file 	gene ID file							##
##	tabular file	gene ID and corrosponding protein IDs				##
##	protein file	protein fasta sequences						##
##											##
##	OUTPUT 	fasta sequences for the corrosponding csv file				##
##											##
## 	eg: python csvGeneID2fasta.py Croaker_DESeq2-results-with-normalized.csv 	##
##		ProteinTable12197_229515.txt GCF_000972845.1_L_crocea_1.0_protein.faa	##
##											##
##########################################################################################

import sys,re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from array import array



def get_sequence(f_file,gid,out_file):
        f = open(f_file, 'r')
        for seq_record in SeqIO.parse(f,"fasta"):
                if ( gid == seq_record.id):
                        rec = SeqRecord(seq_record.seq,
                                id = seq_record.id,
                                description = "")
                        SeqIO.write(rec, out_file, "fasta")
                        break
        f.close()
        return;


def get_csv_geneID(csv_line):
	ln = csv_line.strip().split(",")
	return ln[1].split("\"")[1].split(":")[1]

def get_seqid(g_id,fasta_out):
	gtf = open(sys.argv[2],"r")
	for gtf_line in gtf:
		gtf_line = gtf_line.split('#', 1)[0]
		if(gtf_line != ""):
			gln = gtf_line.strip().split("\t")
			gid = gln[5]
			if(gid == g_id):
				print("GeneID:%s\ttableID:\t%s\tprotein-ID:\t%s" %(g_id, gid, gln[7]))
				get_sequence(sys.argv[3],gln[7],fasta_out)
				break



def main():
	fasta_out = open("fasta_out.fasta", "w")
	f_csv = open(sys.argv[1],"r")
	line = f_csv.next()
	for line in f_csv:
		line = line.split('#', 1)[0]
		if(line != ""):
			gene_id = get_csv_geneID(line)
			##print gene_id			
			get_seqid(gene_id,fasta_out)			
	


main()
