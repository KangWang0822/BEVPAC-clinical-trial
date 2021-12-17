############Data Backup############
/proj/snic2021-23-324/nobackup/private/BEVPAC_WES
sftp kangwang-delivery05046@grus.uppmax.uu.se
get -r UB-2846

############nextflow###############
#data: cd /castor/project/proj/data
#nextflow: cd /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1

#1) Follow the detailed instructions for Uppmax/Bianca (apart from installing Nextflow and nf-core; see step 2): https://github.com/nf-core/configs/blob/master/docs/uppmax.md
#(if you have installed them, you need to undo any changes, and delete any files - to avoid any conflicts)
#2) Load necessary modules: module load bioinfo-tools Nextflow nf-core
#3) Export the following ENV variables, or add them to your .bashrc file:
#4) Run a test run with one paired sample e.g. from WES, using the target BED file, and e.g. the PON from GATK resources. Here is an example:
#module load bioinfo-tools Nextflow nf-core

#!/bin/bash -l
#SBATCH -A sens2019581
#SBATCH -p core 
#SBATCH -n 24
#SBATCH -t 6:00:00
#SBATCH -J sarak_test
module load bioinfo-tools Nextflow nf-core
cd /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"

nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/testdata/test2.tsv" \
--genome GRCh38 \
--germline_resource "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz" \
--germline_resource_index "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi" \
--generate_gvcf \
--pon "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz" \
--pon_index "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi" \
--step 'mapping' \
--tools 'mutect2,strelka,manta,haplotypecaller,cnvkit,snpEff' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/SeqCap_EZ_Exome_v2_hg38_targets.v3.pad.sort.merge.bed" \
--save_bam_mapped \
-resume
