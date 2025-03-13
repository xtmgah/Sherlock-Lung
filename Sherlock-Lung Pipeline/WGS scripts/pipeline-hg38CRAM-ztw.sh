#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with more than one
# set of input fastq files (4 in this example)
# *******************************************
set -e 

# BAM or CRAM file
bamfile=$1
barcode=${bamfile##*/}
#barcode=${barcode%%.*}
barcode=${barcode%%.bam}
barcode=${barcode%%.cram}
foldname=$2
nt=$3
shortReadremove=${4:-false}

ulimit -s 10240

# ******************************************
# 0. Setup
# ******************************************
workdir="/lscratch/$SLURM_JOB_ID/$foldname"
#workdir="/data/Choi_lung/ZTW/AM_hg38_CRAM/$foldname/tmp"
workdir2="$(pwd)/$foldname"
mkdir -p $workdir
mkdir -p $workdir2
logfile=$workdir2/run.log
exec >$logfile 2>&1
cd $workdir

## RG information 
## remove DS long sentense with space. DS can be the last one. 
samtools view -H $bamfile|grep "^@RG"|sed -e 's/\tDS:[^\t]*\([\t\.]\)/\t/' -e  's/[\t]+/\t/g' -e  's/\t$//g'|sed 's/\t/\\t/g' |sed 's/ /_/g'  >RG_info.txt 

## split bam by @RG ## 
samtools split -@ $nt -u norg_file.bam $bamfile

## index all bam ## 
ls ${barcode}_*.bam |rush --dry-run -k -v nt=$nt 'samtools index -@ {nt} {}' |sh
ls -hl

#revise the rg information
#ls ${barcode}_*.bam |sed 's/_\([0-9]*\)\.bam/\t\1/g'|rush --dry-run -k '@RG\tID:LB_{2}\tSM:{1}\tPL:ILLUMINA' >RG_info.txt

## covert to fq file
ls ${barcode}_*.bam|rush --dry-run -k -v nt=$nt 'java -jar $BAZAMPATH/bazam.jar -n {nt} -bam {} |bgzip >{.}.fq.gz && rm -rf {} {}.bai' |sh
ls -hl 

## remove split bam file
rm -rf  ${barcode}_*.bam* 

# remove paired reads with less than 30bp 

if [ "shortReadremove" = true ]; then
  ls *.fq.gz |rush -k -v p="'"  'zcat {}|paste - - - - - - - - | awk -F "\t" -v OFS="\t" {p}{if(length($2)>30 && length($6)>30) print $1"\n"$2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8 |"bgzip >{..}.fastq.gz"; else print $1"\n"$2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8 |"bgzip >{..}_lowq.fastq.gz";}{p} '
  rm -rf *.fq.gz 

else
  for file in *.fq.gz; do
    mv "$file" "${file%.fq.gz}.fastq.gz"
  done

fi

ls -hl 

# Update with the fullpath location of your sample fastq
#fastq_folder=$1
#barcode=$2
#cramoutput=$3

#num_sets=`ls ${fastq_folder}/*_1.fq.gz |wc -l`
# fastq should be named $fastq_prefix$set_number$fastq_suffix_1
# for this case set1_1.fastq.gz, set1_2.fastq.gz, set2_1.fastq.gz, set2_2.fastq.gz, set3_1.fastq.gz, set3_2.fastq.gz
fastq_prefix=$barcode
fastq_suffix=".fastq.gz"
#group_prefix="read_group_name"
#sample="sample_name"
#platform="illumina"
#center="novogene"

# Update with the location of the reference data files
fasta="/data/zhangt8/Ref/Reference/hg38/Homo_sapiens_assembly38.fasta"
dbsnp="/fdb/GATK_resource_bundle/hg38/dbsnp_138.hg38.vcf.gz"
known_Mills_indels="/fdb/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_1000G_indels="/fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"

# Update with the location of the Sentieon software package and license file
release_dir=/data/zhangt8/Ref/Sentieon/sentieon-genomics-202308.03
export SENTIEON_INSTALL_DIR=$release_dir
export SENTIEON_LICENSE=XXXXXXXXXX
export SENTIEON_AUTH_DATA=XXXXXXXXXX

# Other settings
#nt=$(nproc) #number of threads to use in computation, set to number of cores in the server
#nt=16



# ******************************************
# 1. Mapping each set of input fastq with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 

# total fq file
itot=`ls ${barcode}*.fastq.gz|grep -v "lowq.fastq.gz" |wc -l `
itot=$(( itot - 1 ))


i=0
for (( i=0; i<=$itot; i++ )); 
do
  filei=${barcode}"_"$i$fastq_suffix
  rgi=`awk -v line=$i 'NR==(line+1)' RG_info.txt` 
  echo $filei
  echo $rgi
  bam_input="$bam_input -i sorted_set${i}.bam"
 ( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -p -R $rgi -t $nt -K 10000000 $fasta $filei || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $fasta -o sorted_set$i.bam -t $nt --sam2bam -i -
 rm -rf $filei

done

ls -hl 

# ******************************************
# 2. Metrics on the multiple sorted BAM files
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt $bam_input --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
#$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
#$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
#$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
#$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads on the multiple
# sorted BAM files. It is possible
# to mark instead of remove duplicates
# by ommiting the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt $bam_input --algo LocusCollector --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt $bam_input --algo Dedup  --score_info score.txt --metrics dedup_metrics.txt $workdir/deduped.bam 
rm -rf sorted_set*.bam* 

# ******************************************
# 2a. Coverage metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $workdir/deduped.bam --algo CoverageMetrics coverage_metrics

# ******************************************
# 4b. Indel realigner for normal sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $workdir/deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels $workdir/realigned.bam
rm $workdir/deduped.bam* 
# ******************************************
# 5. Base recalibration
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $workdir/realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels $workdir/recal_data.table
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post --algo ReadWriter --cram_write_options compression_level=9,compressor=bzip2 $workdir2/${barcode}.cram
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $workdir2/realigned.bam --read_filter QualCalFilter,table=$workdir2/recal_data.table,indel=false,keep_oq=true,levels=10/20/30/40/50 --algo ReadWriter $workdir2/${barcode}.cram
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $workdir2/realigned.bam --read_filter QualCalFilter,table=$workdir2/recal_data.table,indel=false,levels=10/20/30/40/50 --algo ReadWriter $workdir2/${barcode}.oq.cram
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $workdir2/realigned.bam --read_filter QualCalFilter,table=$workdir2/recal_data.table,indel=true,levels=10/20/30/40/50 --algo ReadWriter $workdir2/${barcode}.ID.cram
## oq generated smallest file size and Mutational calling are almost identical
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $workdir/realigned.bam --read_filter QualCalFilter,table=$workdir/recal_data.table,indel=false,levels=10/20/30/40/50 --algo ReadWriter $workdir2/${barcode}.cram

#$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
#$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o $workdir2/recal_plots.pdf recal.csv

rm $workdir/realigned.bam*
rm $workdir/recal_data.table*

# ******************************************
# 6b. HC Variant caller
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=30 --call_conf=30 output-hc.vcf.gz

ls -lh
### cop file back to workdir ### 
mv RG_info.txt norg_file.bam  $workdir2/

if [ "shortReadremove" = true ]; then
  mv *_lowq.fastq.gz $workdir2/
fi



