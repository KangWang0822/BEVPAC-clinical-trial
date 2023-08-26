#################Re-run Mutect2#######################
cd /proj/sens2019581/nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/Mutect2_germline_added
module load bioinfo-tools 
module load Nextflow/21.04.1 
module load nf-core/1.14 
module load iGenomes/latest
export NXF_OFFLINE='TRUE'
export NXF_HOME="/castor/project/proj/nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/"
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_SINGULARITY_CACHEDIR="/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/"
nextflow run /castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/workflow/main.withgermline.nf \
-profile uppmax \
-with-singularity "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/singularity-images/nf-core-sarek-2.7.1.simg" \
--custom_config_base "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/configs" \
-c "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-data/sarek.custom.config" \
--project sens2019581 \
--input "/proj/sens2019581/nobackup/nf-core2/nf-core-sarek-2.7.1/Sarek-results/Sarek-BEVPAC_WES/results/Preprocessing/TSV/recalibrated.tsv" \
--genome GRCh38 \
--germline_resource "$IGENOMES_DATA/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz" \
--germline_resource_index "$IGENOMES_DATA/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz.tbi" \
--pon "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz" \
--pon_index "/castor/project/proj_nobackup/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi" \
--step 'variant_calling' \
--tools 'mutect2' \
--target_bed "/castor/project/proj_nobackup/nf-core2/nf-core-sarek-2.7.1/BEVPAC-data/Twist_Exome_RefSeq_targets_hg38_100bp_padding.bed" \
-resume
