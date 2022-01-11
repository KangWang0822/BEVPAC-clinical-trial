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
cp -r /castor/project/proj/data/BEVPAC_WES/UB-2845/210820_A00181_0332_AHGCGFDSX2/SampleSheet.csv /proj/nobackup/sens2019581/wharf/kangwang/kangwang-sens2019581
cp -r /proj/nobackup/sens2019581/wharf/kangwang/kangwang-sens2019581/BEVPAC_samplesheet.tsv /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data
cp -r /proj/nobackup/sens2019581/wharf/kangwang/kangwang-sens2019581/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data

#!/bin/bash -l
#SBATCH -A sens2019581
#SBATCH -p core -n 8
#SBATCH -t 24:00:00
#SBATCH -J rm_work
cd /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/
rsync -a --delete work/

###To run nextflow
screen -S sarek
control+A+D
exit 
screen ls
screen -r
screen -X -S 36740.sarek quit



#!/bin/bash -l
#SBATCH -A sens2019581
#SBATCH -p node -n 1 -C mem128GB
#SBATCH -t 240:00:00
#SBATCH -J sarak_test
module load bioinfo-tools Nextflow/21.04.1 nf-core/1.14 iGenomes/latest
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
# Mapping step
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/BEVPAC_samplesheet.tsv" \
--genome GRCh38 \
--step 'mapping' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
--save_bam_mapped \
-resume lethal_engelbart 

-resume [run-name]
#Use the "nextflow log" command to show previous run names.


cd /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES

# Prepare recalibration
module load bioinfo-tools Nextflow/21.04.1 nf-core/1.14 iGenomes/latest
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results/Preprocessing/TSV/duplicates_marked_no_table.tsv" \
--genome GRCh38 \
--step 'prepare_recalibration' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
--save_bam_mapped \
-resume jovial_sanger

# Recalibrate
module load bioinfo-tools Nextflow/21.04.1 nf-core/1.14 iGenomes/latest
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results/Preprocessing/TSV/duplicates_marked.tsv" \
--genome GRCh38 \
--step 'recalibrate' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
--save_bam_mapped \
-resume distraught_poincare

#generate index
cd /castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle
module load tabix
tabix -p vcf 1000g_pon.hg38.vcf.gz

# Variant calling: MuTect2
module load bioinfo-tools Nextflow/21.04.1 nf-core/1.14 iGenomes/latest
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results/Preprocessing/TSV/recalibrated.tsv" \
--genome GRCh38 \
--germline_resource "$IGENOMES_DATA/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz" \
--germline_resource_index "$IGENOMES_DATA/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi" \
--pon "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz" \
--pon_index "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi" \
--step 'variant_calling' \
--tools 'mutect2' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
-resume pensive_bhabha


# Variant calling: Manta (first run Manta, then Strelka ;according to Best Practices)
module load bioinfo-tools Nextflow/21.04.1 nf-core/1.14 iGenomes/latest
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results/Preprocessing/TSV" \
--genome GRCh38 \
--step 'variant_calling' \
--tools 'manta' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
-resume

# Variant calling: Strelka (first run Manta, then Strelka ;according to Best Practices)
module load bioinfo-tools Nextflow/21.04.1 nf-core/1.14 iGenomes/latest
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results/Preprocessing/TSV" \
--genome GRCh38 \
--step 'variant_calling' \
--tools 'strelka' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
-resume

# Variant calling: CNVkit
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.nf \
--outdir "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results" \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results/Preprocessing/TSV" \
--genome GRCh38 \
--step 'variant_calling' \
--tools 'cnvkit' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
-resume





