#!/bin/bash
set -e 
##### Software Variant #######
#module load samtools
#module load GATK
#module load snpEff 
#module load bedtools
#module load annovar 
##### Input ########

mutect2_VCF=$1
tnscope_VCF=$2
strelka_indel_VCF=$3
normal_bam=$4
tumor_bam=$5
PREFIX=$6

fasta="/data/zhangt8/APOBEC/Ref/Homo_sapiens_assembly19.fasta"
upsindel="/data/zhangt8/Ref/ups-indel"
nset=3

DIR=$(pwd)/Result/${PREFIX}

mkdir -p $DIR

cd $DIR

#### extract normal and tumor name ####

normal_name=`samtools view -H ${normal_bam} |grep "@RG" |sed 's/.*SM://'|awk '{print $1}' |uniq`
tumor_name=`samtools view -H  ${tumor_bam} |grep "@RG" |sed 's/.*SM://'|awk '{print $1}' |uniq `

echo "NORMAL $normal_name" >newname.txt
echo "TUMOR $tumor_name" >>newname.txt 

### process orignal vcf call ######
bcftools view -f "PASS,." -v indels -m2 -M2 $mutect2_VCF  >${PREFIX}-tnhaplotyper.vcf 
$upsindel/ups_indel $fasta  ${PREFIX}-tnhaplotyper.vcf ${PREFIX}-tnhaplotyper -hd=true 

bcftools view -f "PASS,." -v indels -m2 -M2  $tnscope_VCF  >${PREFIX}-tnscope.vcf
$upsindel/ups_indel $fasta ${PREFIX}-tnscope.vcf ${PREFIX}-tnscope -hd=true


bcftools reheader -s newname.txt $strelka_indel_VCF |bcftools view -f "PASS,." -v indels -m2 -M2  >${PREFIX}-strelka.vcf
$upsindel/ups_indel $fasta ${PREFIX}-strelka.vcf ${PREFIX}-strelka  -hd=true 

## find overlap 
grep -v "^##" ${PREFIX}-tnscope.uvcf|awk -F "\t" -v OFS="\t" '{print $1,$8,"tnscope"}' >${PREFIX}-indels.list  
grep -v "^##" ${PREFIX}-tnhaplotyper.uvcf|awk -F "\t" -v OFS="\t" '{print $1,$8,"tnhaplotyper"}' >>${PREFIX}-indels.list  
grep -v "^##" ${PREFIX}-strelka.uvcf|awk -F "\t" -v OFS="\t" '{print $1,$8,"strelka"}' >>${PREFIX}-indels.list  

cat ${PREFIX}-indels.list |sort -k1,1 -k2,2n |groupBy -g 1,2 -c 3,3 -o collapse,count_distinct >${PREFIX}-indels.list2
mv ${PREFIX}-indels.list2 ${PREFIX}-indels.list
cat ${PREFIX}-indels.list |cut -f 3 |sort|uniq -c >${PREFIX}-indels.list.overlap


## filter orignal vcf 
cat ${PREFIX}-indels.list |awk -F "\t" -v OFS="\t" -v nset=$nset '$4>=nset' |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1"@"$2]=1;next;} ($1"@"$8 in a) {print $0}' - ${PREFIX}-tnhaplotyper.uvcf  |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1"@"$2"@"$4"@"$5]=1;next;}/^#/ || ($1"@"$2"@"$4"@"$5) in a{print $0}' - ${PREFIX}-tnhaplotyper.vcf >${PREFIX}-tnhaplotyper-filter.vcf 
cat ${PREFIX}-indels.list |awk -F "\t" -v OFS="\t" -v nset=$nset '$4>=nset' |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1"@"$2]=1;next;} ($1"@"$8 in a) {print $0}' - ${PREFIX}-tnscope.uvcf  |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1"@"$2"@"$4"@"$5]=1;next;}/^#/ || ($1"@"$2"@"$4"@"$5) in a{print $0}' - ${PREFIX}-tnscope.vcf >${PREFIX}-tnscope-filter.vcf 
cat ${PREFIX}-indels.list |awk -F "\t" -v OFS="\t" -v nset=$nset '$4>=nset' |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1"@"$2]=1;next;} ($1"@"$8 in a) {print $0}' - ${PREFIX}-strelka.uvcf  |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1"@"$2"@"$4"@"$5]=1;next;}/^#/ || ($1"@"$2"@"$4"@"$5) in a{print $0}' - ${PREFIX}-strelka.vcf >${PREFIX}-strelka-filter.vcf 


## combine and normalize variants
java -Xmx4g -jar  $GATK_JAR  -T CombineVariants -R  $fasta  -genotypeMergeOptions PRIORITIZE --variant:mutect2 ${PREFIX}-tnhaplotyper-filter.vcf --variant:tnscope ${PREFIX}-tnscope-filter.vcf --variant:strelka ${PREFIX}-strelka-filter.vcf -o combine_indel.vcf -priority mutect2,tnscope,strelka  --minimumN $nset --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED

bcftools norm -f $fasta -o combine_indel_lt.vcf combine_indel.vcf 

### filtered by QC ###

fname1=`bcftools view -h combine_indel_lt.vcf |tail -1|cut -f 10` 
fname2=`bcftools view -h combine_indel_lt.vcf |tail -1|cut -f 11`

if [ "$fname1" != "$tumor_name" ]
then
	echo "Tumor ID in vcf file is not the first"
	java -jar $SNPSIFT_JAR extractFields combine_indel_lt.vcf CHROM POS REF ALT GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP GEN[0].AD[0] GEN[0].AD[1] GEN[0].DP >${PREFIX}.VAF
	#exit
else 
	echo "Tumor ID in vcf file is the first"
	java -jar $SNPSIFT_JAR extractFields  combine_indel_lt.vcf CHROM POS REF ALT GEN[0].AD[0] GEN[0].AD[1] GEN[0].DP  GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP >${PREFIX}.VAF
fi 

cat ${PREFIX}.VAF|awk -F "\t" -v OFS="\t" 'NR>1{if(($5==0 && $6==0)||($8==0 && $9==0)) print $1,$2,$3,$4,"0","0"; else {TAF=$6/($5+$6);NAF=$9/($8+$9); if(NAF>0.02||$6<3||($5+$6)<8||($8+$9)<6) print $1,$2,$3,$4,TAF,NAF}}' >${PREFIX}.VAF.bad
#TAF<0.04

if [ -s ${PREFIX}.VAF.bad ]
then
		cat combine_indel_lt.vcf | awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1":"$2":"$3":"$4]=1;next}{if(!($1":"$2":"$4":"$5 in a)) print $0 }' ${PREFIX}.VAF.bad  -  >${PREFIX}-final-indels.vcf
else
        echo "${PREFIX}.VAF.bad is empty!!!"
        cp combine_indel_lt.vcf  ${PREFIX}-final-indels.vcf 
fi

mv combine_indel_lt.vcf ${PREFIX}-inital-indels.vcf 

## annotation format ##
convert2annovar.pl -format vcf4 ${PREFIX}-final-indels.vcf -outfile final-indels -allsample

mv final-indels.${tumor_name}.avinput ${PREFIX}-final-indels.avinput 
## clear folder ####

ls ./|grep -v "final-indels.vcf" |grep -v "final-indels.avinput"|grep -v "overlap" |grep -v "inital-indels.vcf"|awk '{print "rm -rf "$1}'|sh 



#convert2annovar.pl -format vcf4 $VCF -outfile ${Sample}.avinput

table_annovar.pl ${PREFIX}-final-indels.avinput $ANNOVAR_DATA/hg19 \
  --tempdir /lscratch/$SLURM_JOB_ID \
  --thread 8 \
  --buildver hg19 \
  --outfile ${PREFIX}-final-indels \
  --remove \
  --protocol refGene,ensGene,avsift,ljb26_all,dbnsfp30a,cg46,dbscsnv11,cosmic64,cosmic70,exac03,exac03nontcga,1000g2015aug_all,1000g2012apr_all,snp138,avsnp147,clinvar_20160302,gnomad_exome,gnomad_genome \
  --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
  --nastring ''
