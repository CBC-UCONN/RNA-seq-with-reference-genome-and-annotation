#!/bin/bash

# run all tutorial scripts using slurm dependencies

cd 01_raw_data
jid1=$(sbatch --parsable 01_fasterq_dump.sh)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 02_get_genome.sh)

cd ../02_quality_control/
jid3=$(sbatch --parsable --dependency=afterok:$jid2 01_fastp_trim.sh)
jid4=$(sbatch --parsable --dependency=afterok:$jid3 02_fastqc.sh)
jid5=$(sbatch --parsable --dependency=afterok:$jid4 03_multiqc.sh)

cd ../03_index/
jid6=$(sbatch --parsable --dependency=afterok:$jid5 hisat2_index.sh)

cd ../04_align/
jid7=$(sbatch --parsable --dependency=afterok:$jid6 align.sh)

cd ../05_align_qc/
jid8=$(sbatch --parsable --dependency=afterok:$jid7 01_samtools_stats.sh)
jid9=$(sbatch --parsable --dependency=afterok:$jid8 02_qualimap.sh)
jid10=$(sbatch --parsable --dependency=afterok:$jid9 03_multiqc.sh)

cd ../06_count/
jid11=$(sbatch --parsable --dependency=afterok:$jid10 htseq_count.sh)

cd ../07_stringtie/
jid12=$(sbatch --parsable --dependency=afterok:$jid11 01_stringtie_assemble.sh)
jid13=$(sbatch --parsable --dependency=afterok:$jid12 02_stringie_merge.sh)
jid14=$(sbatch --parsable --dependency=afterok:$jid13 03_gtf_annotate.sh)

cd ../08_kallisto
jid15=$(sbatch --parsable --dependency=afterok:$jid14 01_kallisto_stringtie.sh)
jid16=$(sbatch --parsable --dependency=afterok:$jid15 02_kallisto_ensembl.sh)
