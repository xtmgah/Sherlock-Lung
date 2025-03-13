#!/bin/sh
# *******************************************
# Script to perform TN seq variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************
set -e 
##modify by ZTW##
tumor_sample=$1
tumor_bam=$2
normal_sample=$3
normal_bam=$4
lastname=$5
rgfile=$6

ulimit -s 10240
#BAM_folder=$3
#BAM_folder="/data/zhangt8/AM/AM_China_WGS/Sentieon/BAMs"

# Update with the location of the reference data files
#fasta="/data/zhangt8/APOBEC/Ref/Homo_sapiens_assembly19.fasta"
fasta="/data/zhangt8/Ref/Reference/hg38/Homo_sapiens_assembly38.fasta"
dbsnp="/fdb/GATK_resource_bundle/hg38bundle/dbsnp_146.hg38.vcf.gz"
known_Mills_indels="/fdb/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_1000G_indels="/fdb/GATK_resource_bundle/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf.gz" 

# Update with the location of the Sentieon software package and license file
release_dir=/data/zhangt8/Ref/Sentieon/sentieon-genomics-202308.03
export SENTIEON_LICENSE=XXXXXXXX
export SENTIEON_AUTH_DATA=XXXXXXXX

# Other settings
nt=16 #number of threads to use in computation
workdir="/lscratch/$SLURM_JOB_ID/$lastname"
workdir2=`pwd`"/TNseq/$lastname"  #Determine where the output files will be stored
# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
mkdir -p $workdir2
logfile=$workdir2/run.log
exec >$logfile 2>&1
cd $workdir

## extract tumor and normal sample name 
#tumor_sample=`samtools view -H $tumor_bam |grep "^@RG" |sed 's/.*\tSM://g'|awk '{print $1}' |uniq |head -1` 
#normal_sample=`samtools view -H $normal_bam |grep "^@RG" |sed 's/.*\tSM://g'|awk '{print $1}' |uniq|head -1 `

## check the IDs between tumro and normal
#samtools view -H $normal_bam |grep "^@RG" |sed 's/.*\tID://g'|awk '{print $1}' |sort|uniq >RG.txt
#samtools view -H $tumor_bam |grep "^@RG" |sed 's/.*\tID://g'|awk '{print $1}' |sort|uniq >>RG.txt
#dupID=`cat RG.txt |sort|uniq -d|wc -l`
#if [[ $dupID != 0 ]]; then
#  samtools view -h $normal_bam | sed -e 's/\(^@RG\tID:\)\([^\s]*\)/\1N_\2/' -e 's/\t\(RG:Z:\)/\t\1N_/g' | samtools view -b -@ $nt -o normal.bam -
#  samtools index -@ $nt normal.bam
#  normal_bam="normal.bam"
#fi


# ******************************************
# 6. Corealignment of tumor and normal
# ******************************************

# check the Read group information
( samtools view -H $tumor_bam ; samtools view -H $normal_bam ) | grep "^@RG" | sed -n 's/^@RG.*ID:\([^[:space:]]*\).*SM:\([^[:space:]]*\).*/\2\t\1/p' >all_rgID.txt
echo "### read group information ###" 
cat all_rgID.txt

rgrepeat=$( cat all_rgID.txt |cut -f 2|sort|uniq -d )

addparm=" "

if [ -z "$rgrepeat" ]; then
  
  echo "Great! no duplicate RG ID found!"

else

  if [ -z "$rgfile" ]; then
    echo "Failed, A rg file need to provide"
    exit 0
  else
    addparm=$(cat $rgfile | sed 's/.*ID://g' | sed 's|\\.*||g' | paste - $rgfile | sed 's|@RG\\tID:|ID:|' | rush --dry-run -k " --replace_rg {1}='{2}'" | sed 's/ID:/ID:N_/' | paste -s -d " ")
  fi
  
fi


echo "###### $addparm #####" 

eval "$release_dir/bin/sentieon driver -r $fasta  -t $nt -i  $tumor_bam $addparm -i $normal_bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tn_corealigned.bam" 

rm all_rgID.txt

# ******************************************
# 7. Somatic Variant Calling
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tn_corealigned.bam --algo TNsnv --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp --call_stats_out ${lastname}-call.stats ${lastname}-tnsnv.vcf.gz
cp $workdir/${lastname}-tnsnv.vcf.gz* $workdir2 

$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tn_corealigned.bam --algo TNhaplotyper --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp ${lastname}-tnhaplotyper.vcf.gz
cp $workdir/${lastname}-tnhaplotyper.vcf.gz* $workdir2 

# ******************************************
# 8. Somatic and Structural variant calling
# ******************************************
$release_dir/bin/sentieon driver -r $fasta  -t $nt -i tn_corealigned.bam --algo TNscope --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp ${lastname}-tnscope.vcf.gz
cp $workdir/${lastname}-tnscope.vcf.gz*  $workdir2 

# ******************************************
# 8a. HC Variant caller for GVCF
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $tumor_bam  --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 --emit_mode gvcf ${tumor_sample}_hc.g.vcf.gz
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $normal_bam  --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=30 --emit_mode gvcf ${normal_sample}_hc.g.vcf.gz

cp $workdir/${tumor_sample}_hc.g.vcf.gz*  $workdir/${normal_sample}_hc.g.vcf.gz* $workdir2 


### calling with DNAscope ### 
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $tumor_bam  --algo DNAscope -d $dbsnp ${tumor_sample}-dnascope.vcf.gz
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $normal_bam  --algo DNAscope -d $dbsnp ${normal_sample}-dnascope.vcf.gz

cp $workdir/${tumor_sample}-dnascope.vcf.gz* $workdir/${normal_sample}-dnascope.vcf.gz* $workdir2 
### Structural variantion calling by DNAscope ##
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $tumor_bam  --algo DNAscope --var_type bnd -d $dbsnp ${tumor_sample}-tmpdna.vcf.gz
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $normal_bam  --algo DNAscope --var_type bnd -d $dbsnp ${normal_sample}-tmpdna.vcf.gz


$release_dir/bin/sentieon driver -r $fasta -t $nt -i $tumor_bam --algo SVSolver -v ${tumor_sample}-tmpdna.vcf.gz ${tumor_sample}-dnascopeSV.vcf.gz 
$release_dir/bin/sentieon driver -r $fasta -t $nt -i $normal_bam --algo SVSolver -v ${normal_sample}-tmpdna.vcf.gz ${normal_sample}-dnascopeSV.vcf.gz 

cp $workdir/${tumor_sample}-dnascopeSV.vcf.gz* $workdir/${normal_sample}-dnascopeSV.vcf.gz* $workdir2

rm tn_corealigned.bam*
