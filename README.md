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
3. Quality control using sickle
4. Aligning Reads to a Genome using hisat2
5. Generating Total Read Counts from Alignment using htseq-count
6. Pairwise differential expression with counts in R with DESeq2
	1. Common plots for differential expression analysis
	2. Using DESeq2
7. EnTAP: Functional Annotation for Genomes
8. Integrating the DE Results with the Annotation Results  


### 1. Overview  

Liver mRNA profiles large yellow croaker (Larimichthys crocea) species are sampled during various conditions namely, control group (LB2A), thermal stress group (LC2A), cold stress group (LA2A) and 21-day fasting group (LF1A) were generated by RNA-seq, using Illumina HiSeq 2000.  

We will only use the control group (LB2A) and the thermal stress group (LC2A) 

The workflow may be cloned into the appropriate directory using the terminal command:  
```bash
git clone https://github.com/golden75/Model_Marine_RNA_Seq_and_functionalAnnotation.git 
```  

### 2. Accessing the Data using SRA-Toolkit    


We will be downloading our data from the sequence-read-archives (SRA), a comprehensive collection of sequenced genetic data submitted to the NCBI by experimenters. The beauty of the SRA is the ease with which genetic data becomes accessible to any scientist with an internet connection, available for download in a variety of formats. Each run in the SRA has a unique identifier. The run may be downloaded using a module of software called the "sratoolkit" and its unique identifier. There are a variety of commands in the sratoolkit, which I invite you to investigate for yourself at https://www.ncbi.nlm.nih.gov/books/NBK158900/.  

The data may be accessed at the following web page: https://www.ncbi.nlm.nih.gov/bioproject/28084  
LB2A : SRR1964642, SRR1964643  
LC2A : SRR1964644, SRR1964645  

We can download the SRA data using SRA-Toolkit using the following command:  
```bash
#!/bin/bash
#SBATCH --job-name=fastq_dump_xanadu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=15G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

module load sratoolkit/2.8.1  
fastq-dump SRR1964642  
fastq-dump SRR1964643
fastq-dump SRR1964644  
fastq-dump SRR1964645  
```

Once the download is complete, we will rename the according to the samples for easy identification.  
```bash
mv SRR1964642.fastq LB2A_SRR1964642.fastq
mv SRR1964643.fastq LB2A_SRR1964643.fastq
mv SRR1964644.fastq LC2A_SRR1964644.fastq
mv SRR1964645.fastq LC2A_SRR1964645.fastq
```  

The full script for slurm scheduler can be found in the raw_data folder by the name [fastq_dump_xanadu.sh](/raw_data/fastqc_dump_xanadu.sh).  

It is advised that you familiarize yourself with the arguments for the Slurm scheduler. While it may seem as though running your commands locally will be more efficient due to the hassle of not initializing and writing scripts, do not fall for that trap! The capacity of the Slurm scheduler far exceeds the quickness of entering the commands locally. While the rest of this tutorial will not include the process of initializing and writing the Slurm arguments in a script in its coding, know that the Xanadu scripts in the cloned directory do contain the Slurm arguments. However, before running any cloned Xanadu script, you must edit and enter your appropriate email address!  

Once the download and the renaming is done the folder structure will look like:  
```bash
raw_data/
|-- LB2A_SRR1964642.fastq
|-- LB2A_SRR1964643.fastq
|-- LC2A_SRR1964644.fastq
`-- LC2A_SRR1964645.fastq
```   

Lets have a look at athe contents of one of the fastq-files:  
```bash
head -n 12 LB2A_SRR1964642.fastq

@SRR1964642.1 FCC355RACXX:2:1101:1476:2162 length=90
CAACATCTCAGTAGAAGGCGGCGCCTTCACCTTCGACGTGGGGAATCGCTTCAACCTCACGGGGGCTTTCCTCTACACGTCCTGTCCGGA
+SRR1964642.1 FCC355RACXX:2:1101:1476:2162 length=90
?@@D?DDBFHHFFGIFBBAFG:DGHDFHGHIIIIC=D<:?BBCCCCCBB@BBCCCB?CCBB<@BCCCAACCCCC>>@?@88?BCACCBB>
@SRR1964642.2 FCC355RACXX:2:1101:1641:2127 length=90
NGCCTGTAAAATCAAGGCATCCCCTCTCTTCATGCACCTCCTGAAATAAAAGGGCCTGAATAATGTCGTACAGAAGACTGCGGCACAGAC
+SRR1964642.2 FCC355RACXX:2:1101:1641:2127 length=90
#1=DDFFFHHHHGJJJJJIIIJIJGIIJJJIJIJJGIJIJJJJIJJJJJJIJJJIJJJJJJJGIIHIGGHHHHHFFFFFDEDBDBDDDDD
@SRR1964642.3 FCC355RACXX:2:1101:1505:2188 length=90
GGACAACGCCTGGACTCTGGTTGGTATTGTCTCCTGGGGAAGCAGCCGTTGCTCCACCTCCACTCCTGGTGTCTATGCCCGTGTCACCGA
+SRR1964642.3 FCC355RACXX:2:1101:1505:2188 length=90
CCCFFFFFHHFFHJJJIIIJHHJJHHJJIJIIIJEHJIJDIJJIIJJIGIIIIJGHHHHFFFFFEEEEECDDDDEDEDDDDDDDADDDDD
```  

We see that for our first three runs we have information about the sampled read including its length followed by the nucleotide read and then a "+" sign. The "+" sign marks the beginning of the corresponding scores for each nucleotide read for the nucleotide sequence preceding the "+" sign.  


## 3. Quality Control Using Sickle  


