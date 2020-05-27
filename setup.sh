#!/bin/bash

# run all tutorial scripts using slurm dependencies

cd raw_data
jid1=$(sbatch --parsable fastq_dump_xanadu.sh)

cd ../quality_control/
jid2=$(sbatch --parsable --dependency=afterok:$jid1 fastq_trimming.sh)

cd ../fastqc/
jid3=$(sbatch --parsable --dependency=afterok:$jid2 fastqc_raw.sh)
jid4=$(sbatch --parsable --dependency=afterok:$jid3 fastqc_trimmed.sh)
jid5=$(sbatch --parsable --dependency=afterok:$jid4 multiqc_raw.sh)
jid6=$(sbatch --parsable --dependency=afterok:$jid5 multiqc_trimmed.sh)

cd ../index/
jid7=$(sbatch --parsable --dependency=afterok:$jid6 hisat2_index.sh)

cd ../align/
jid8=$(sbatch --parsable --dependency=afterok:$jid7 align.sh)

cd ../count/
jid9=$(sbatch --parsable --dependency=afterok:$jid8 htseq_count.sh)
