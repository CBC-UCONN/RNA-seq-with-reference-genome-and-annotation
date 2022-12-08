#!/bin/bash 
#SBATCH --job-name=samstats
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

##################################
# calculate stats on alignments
##################################

# calculate alignment statistics for each bam file using samtools

# load software--------------------------------------------------------------------------
module load samtools/1.16.1
module load parallel/20180122

# input, output directories--------------------------------------------------------------

INDIR=../04_align/alignments
OUTDIR=samtools_stats
mkdir -p $OUTDIR

ACCLIST=../01_raw_data/accessionlist.txt

# samtools bam statistics----------------------------------------------------------------

cat $ACCLIST | parallel -j 5 \
	"samtools stats $INDIR/{}.bam >$OUTDIR/{}.stats"

# put the basic stats all in one file.---------------------------------------------------

# bash array containing all stats file names
FILES=($(find $OUTDIR -name "*.stats" | sort))

# initialize big table for all samples with row names
grep "^SN" ${FILES[0]} | cut -f 2 > $OUTDIR/SN.txt

# loop over each stats file, add data column to main table
for file in ${FILES[@]}
do paste $OUTDIR/SN.txt <(grep ^SN $file | cut -f 3) > $OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt
	echo $file
done

# add a header with sample names
cat \
<(echo ${FILES[@]} | sed 's,samtools_stats/,,g' | sed 's/.stats//g' | sed 's/ /\t/g') \
$OUTDIR/SN.txt \
>$OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt