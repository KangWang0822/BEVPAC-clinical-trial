#BEVPAC clinical trial

#Data deliver for NGI
#Data location
cd /proj/sens2019581/delivery04273/INBOX/P18362/

#传输：https://www.uppmax.uu.se/support/user-guides/basic-sftp-commands/
#Uploaded 
sftp -q kangwang-sens2019581@bianca-sftp.uppmax.uu.se
cd kangwang-sens2019581
lls
put hg38_chr.map  #上传
##登录terminal
cp hg38_chr.map /home/kangwang/tools/CrosscheckFingerprints
cp -r /home/kangwang
#Download
cd /proj/sens2019581/nobackup/wharf/kangwang/kangwang-sens2019581 #All files must be transferred through the wharf area of Bianca (download) 复制文件到该目录下
ls


#SBATCH
squeue -u username
sbatch name.slurm

#creat a file and write bash script
cat >name.sh
content.....
Ctrl+D  #ending edition
Vim name.sh # editing modify

################################
################################
################################
#!/bin/bash -l

#SBATCH -A sens2019581
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -J CrosscheckFingerprints 

biopsy1="133"
biopsy2="165"

IN1=$(ls /proj/sens2019581/delivery04273/INBOX/P18362/01-RNA-Results/markDuplicates/P18362_${biopsy1}Aligned.sortedByCoord.out.markDups.bam)
IN2=$(ls /proj/sens2019581/delivery04273/INBOX/P18362/01-RNA-Results/markDuplicates/P18362_${biopsy2}Aligned.sortedByCoord.out.markDups.bam)

echo "Loading necessary modules..."
module load bioinfo-tools picard/2.23.4

echo "Running CrosscheckFingerprints for the paired samples ${biopsy1} and ${biopsy2}..."
java -jar $PICARD_ROOT/picard.jar CrosscheckFingerprints \
           INPUT=${IN1} \
           INPUT=${IN2} \
           HAPLOTYPE_MAP=/home/kangwang/tools/CrosscheckFingerprints/hg38_chr.map \
           LOD_THRESHOLD=-5 \
	   NUM_THREADS=16 \
           EXPECT_ALL_GROUPS_TO_MATCH=true \
           OUTPUT=/home/kangwang/tools/CrosscheckFingerprints/sample.crosscheck_metrics_${biopsy1}_vs${biopsy2}
echo "Done..."



export NXF_OFFLINE='/proj/sens2019581'
