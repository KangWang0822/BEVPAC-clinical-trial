module load bioinfo-tools
module load Nextflow

touch samplesheet.csv
vim samplesheet.csv
replicate,fastq_1,fastq_2,strandedness
1,/proj/sens2019581/delivery04273/INBOX/P18362/P18362_131/02-FASTQ/210122_A00187_0419_AHNKKTDSXY/P18362_131_S31_L002_R1_001.fastq.gz,/proj/sens2019581/delivery04273/INBOX/P18362/P18362_131/02-FASTQ/210122_A00187_0419_AHNKKTDSXY/P18362_131_S31_L002_R2_001.fastq.gz,unstranded
1,/proj/sens2019581/delivery04273/INBOX/P18362/P18362_109/02-FASTQ/210122_A00187_0419_AHNKKTDSXY/P18362_109_S20_L002_R1_001.fastq.gz,/proj/sens2019581/delivery04273/INBOX/P18362/P18362_109/02-FASTQ/210122_A00187_0419_AHNKKTDSXY/P18362_109_S20_L002_R2_001.fastq.gz,unstranded
按住shift
zz 保存退出

touch genomefastas.csv
vim genomefastas.csv
