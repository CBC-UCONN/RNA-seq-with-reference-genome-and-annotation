# Model Marine RNASeq and Functional Annotation
# RNA-Seq: Reference Genome, Differential Expression, and Functional Annotation

This repository is a usable, publicly available tutorial for analyzing differential expression data and creating topological gene networks. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.

Contents   

	1. Overview and programs install
	2. Accessing the data using sra-toolkit
	3. Quality control using sickle
	4. Aligning reads to a genome using hisat2 
	5. Generating total read counts from alignment using htseq-count 
	6. Pairwise differential expression with counts in R with DESeq2
		1. Common plots for differential expression analysis
		2. Using DESeq2
	7. EnTAP: Functional Annotation for Genomes
	8. Integrating the DE Results with the Annotation Results  

Contents
1. Overview
2. Accessing the Data using SRA-Toolkit  



