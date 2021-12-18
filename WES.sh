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
#SBATCH -p node 
#SBATCH -n 1
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

################################
###########Test for BEV#########
################################
#upload phenotype#
cd /Users/kangwang/KI/Projects/3.BEVPAC
sftp kangwang-sens2019581@bianca-sftp.uppmax.uu.se:kangwang-sens2019581
##
cp -r /proj/nobackup/sens2019581/wharf/kangwang/kangwang-sens2019581/BEVPAC_test.tsv /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data
cp -r /proj/nobackup/sens2019581/wharf/kangwang/kangwang-sens2019581/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data

#!/bin/bash -l
#SBATCH -A sens2019581
#SBATCH -p core 
#SBATCH -n 12
#SBATCH -t 24:00:00
#SBATCH -J sarak_test
module load bioinfo-tools Nextflow/21.04.1 nf-core/1.14
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
# Mapping step
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/BEVPAC_test.tsv" \
--genome GRCh38 \
--step 'mapping' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
--save_bam_mapped \
-resume


# Prepare recalibration

nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/home/emmasi/results/Sarek-mapping/results/Preprocessing/TSV/duplicates_marked_no_table.tsv" \
--genome GRCh38 \
--step 'prepare_recalibration' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/SeqCap_EZ_Exome_v2_hg38_targets.v3.pad.sort.merge.bed" \
--save_bam_mapped \
-resume


# Recalibrate

nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-PROMIX_WES/results/Preprocessing/TSV/duplicates_marked.tsv" \
--genome GRCh38 \
--step 'recalibrate' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/SeqCap_EZ_Exome_v2_hg38_targets.v3.pad.sort.merge.bed" \
--save_bam_mapped \
-resume


# Variant calling: MuTect2

nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-PROMIX_WES/results/Preprocessing/TSV/recalibrated.tsv" \
--genome GRCh38 \
--germline_resource "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz" \
--germline_resource_index "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi" \
--pon "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz" \
--pon_index "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi" \
--step 'variant_calling' \
--tools 'mutect2' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/SeqCap_EZ_Exome_v2_hg38_targets.v3.pad.sort.merge.bed" \
-resume








