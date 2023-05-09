#!/bin/bash
set -e 
##### Input ########

mutect_VCF=$1
mutect2_VCF=$2
tnscope_VCF=$3
strelka_snv_VCF=$4
strelka_indel_VCF=$5
normal_bam=$6
tumor_bam=$7
PREFIX=$8

PREFIX0=$PREFIX

fasta="/data/zhangt8/Ref/Reference/hg38/Homo_sapiens_assembly38.fasta"
source /data/zhangt8/Ref/Reference/hg38/hg38_env.sh

DIR=$(pwd)/Result/${PREFIX}
MAFOUTPUT=${DIR}/$PREFIX}.maf
SGE_LOG=${DIR}/${PREFIX}.log

mkdir -p $DIR
cd $DIR

#### extract normal and tumor name ####
normal_name=`samtools view -T $fasta -H ${normal_bam} |grep "^@RG" |sed 's/.*SM://'|awk '{print $1}' |uniq`
tumor_name=`samtools view -T $fasta -H  ${tumor_bam} |grep "^@RG" |sed 's/.*SM://'|awk '{print $1}' |uniq `

echo "NORMAL $normal_name" >newname.txt
echo "TUMOR $tumor_name" >>newname.txt 

### process snp call ######
PREFIX=${PREFIX0}.snp

bcftools view -i 'SVTYPE="."' $tnscope_VCF -O z -o tnscope_VCF.vcf.gz 
bcftools index -t tnscope_VCF.vcf.gz 

bcftools reheader -s newname.txt $strelka_snv_VCF |bcftools view -O z >strelka.somatic.snvs.vcf.gz
bcftools index -t strelka.somatic.snvs.vcf.gz 
GATK CombineVariants -R $fasta -genotypeMergeOptions PRIORITIZE --variant:mutect2 $mutect2_VCF --variant:mutect $mutect_VCF  --variant:tnscope tnscope_VCF.vcf.gz --variant:strelka strelka.somatic.snvs.vcf.gz -o combine_snv.vcf -priority mutect2,mutect,tnscope,strelka --minimumN 3 --filteredrecordsmergetype KEEP_IF_ALL_UNFILTERED
bcftools view  -m2 -M2 -i ' set="Intersection" | set="mutect2-mutect-tnscope" |set="mutect2-mutect-strelka" | set="mutect2-tnscope-strelka" |set="mutect-tnscope-strelka"  |set="filterInmutect2-mutect-tnscope-strelka" |set="mutect2-filterInmutect-tnscope-strelka"|set="mutect2-mutect-filterIntnscope-strelka"|set="mutect2-mutect-tnscope-filterInstrelka" ' combine_snv.vcf  >combine_snv_filter.vcf 

python /data/zhangt8/mWGS/Squamous_Cell_Carcinoma/SNV/add_to_vcf_header.py combine_snv_filter.vcf '##INFO=<ID=OVERLAP,Number=0,Type=Flag,Description="Somatic indel possibly overlaps a second indel.">' >combine_snv_filter.vcf2
mv combine_snv_filter.vcf2 combine_snv_filter.vcf

#java -cp $GATK_JAR org.broadinstitute.gatk.tools.CatVariants  -R $fasta -V combine_snv_filter.vcf  -V combine_indel_filter.vcf -assumeSorted -out combine_all_filter.vcf 
mv combine_snv_filter.vcf combine_all_filter.vcf  

bcftools sort  -T /lscratch/$SLURM_JOB_ID combine_all_filter.vcf |bcftools view -s ${tumor_name},${normal_name} -o ${PREFIX}.vcf.gz -O z 
bcftools index -f ${PREFIX}.vcf.gz

rm -rf strelka.somatic.snvs.vcf.gz* combine_snv.vcf* combine_snv_filter.vcf* strelka.somatic.indels.vcf.gz* combine_indel.vcf* combine_indel_filter.vcf* combine_all_filter.vcf* tnscope_VCF.vcf.gz*


### process indel call ###
PREFIX=${PREFIX0}.indel
nset_indel=3

bcftools view -f "PASS,." -v indels -m2 -M2 $mutect2_VCF|awk -F "\t" -v OFS="\t" '/^#/ || $5!~/[^ATCGNLT]/' |bcftools norm -f $fasta -o ${PREFIX}-tnhaplotyper.vcf 
bcftools view -f "PASS,." -v indels -m2 -M2  $tnscope_VCF|awk -F "\t" -v OFS="\t" '/^#/ || $5!~/[^ATCGNLT]/' |bcftools norm -f $fasta -o ${PREFIX}-tnscope.vcf
bcftools reheader -s newname.txt $strelka_indel_VCF |bcftools view -f "PASS,." -v indels -m2 -M2|awk -F "\t" -v OFS="\t" '/^#/ || $5!~/[^ATCGNLT]/'  |bcftools norm -f $fasta -o ${PREFIX}-strelka.vcf

GATK CombineVariants -R  $fasta  -genotypeMergeOptions PRIORITIZE --variant:mutect2 ${PREFIX}-tnhaplotyper.vcf --variant:tnscope ${PREFIX}-tnscope.vcf --variant:strelka ${PREFIX}-strelka.vcf -o combine_indel.vcf -priority mutect2,tnscope,strelka  --minimumN $nset_indel --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED
bcftools norm -f $fasta combine_indel.vcf |bcftools view -m2 -M2 -s ${tumor_name},${normal_name} |bcftools sort  -T /lscratch/$SLURM_JOB_ID -O z -o ${PREFIX}.vcf.gz
bcftools index -f ${PREFIX}.vcf.gz

rm -rf ${PREFIX}-tnhaplotyper.vcf* ${PREFIX}-tnscope.vcf* ${PREFIX}-strelka.vcf* combine_indel.vcf* 

### cobmined snp and indel calling ### 
PREFIX=$PREFIX0
bcftools concat -a -O z -o ${PREFIX}.vcf.gz ${PREFIX}.snp.vcf.gz ${PREFIX}.indel.vcf.gz  

mv ${PREFIX}.snp.vcf.gz ${PREFIX}.snp.raw.vcf.gz
mv ${PREFIX}.snp.vcf.gz.csi ${PREFIX}.snp.raw.vcf.gz.csi
mv ${PREFIX}.indel.vcf.gz ${PREFIX}.indel.raw.vcf.gz
mv ${PREFIX}.indel.vcf.gz.csi ${PREFIX}.indel.raw.vcf.gz.csi

### filtered by PASS ###

fname1=`bcftools view -h ${PREFIX}.vcf.gz |tail -1|cut -f 10` 
fname2=`bcftools view -h ${PREFIX}.vcf.gz |tail -1|cut -f 11`

if [ "$fname1" != "$tumor_name" ]
then
	echo "Tumor ID in vcf file is not the first"
	java -jar $SNPSIFT_JAR extractFields  ${PREFIX}.vcf.gz CHROM POS REF ALT GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP GEN[0].AD[0] GEN[0].AD[1] GEN[0].DP >${PREFIX}.VAF
else 
	echo "Tumor ID in vcf file is the first"
	java -jar $SNPSIFT_JAR extractFields  ${PREFIX}.vcf.gz CHROM POS REF ALT GEN[0].AD[0] GEN[0].AD[1] GEN[0].DP  GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP >${PREFIX}.VAF
fi 

### #### filtering steps  ##### 
#cat ${PREFIX}.VAF|awk -F "\t" -v OFS="\t" 'NR>1{if(($5==0 && $6==0)||($8==0 && $9==0)) print $1,$2,$3,$4,"0","0"; else {TAF=$6/($5+$6);NAF=$9/($8+$9); if(NAF>0.02||$6<5||($5+$6)<12||($8+$9)<6) print $1,$2,$3,$4,TAF,NAF}}' >${PREFIX}.VAF.bad
## TAF<0.07| 
#cat ${PREFIX}.VAF|awk -F "\t" -v OFS="\t" 'BEGIN{a1=0;a2=0;a3=0;a4=0;a5=0}NR>1{if(($5==0 && $6==0)||($8==0 && $9==0)) {a1++;print $1,$2,$3,$4,"0","0";} else {TAF=$6/($5+$6);NAF=$9/($8+$9); if(NAF>0.02) a2++; if($6<4) a3++; if(($5+$6)<12) a4++; if(($8+$9)<6) a5++; if(NAF>0.02||$6<5||($5+$6)<12||($8+$9)<6) print $1,$2,$3","$4}}END{print "0_coverage:\t"a1"\nNAF>0.02:\t"a2"\nAD_tumor_var<4:\t"a3"\nAD_tumor:\t"a4"\nAD_normal<6:\t"a5 >"filter_log.txt"}'|bgzip -c  >${PREFIX}.filter1.txt.gz 
#tabix -s1 -b2 -e2 ${PREFIX}.filter1.txt.gz

MNAF=2
MRDT=12
MRDN=6

cat ${PREFIX}.VAF | awk -F "\t" -v OFS="\t" -v f2=$MNAF -v f4=$MRDT -v f5=$MRDN 'BEGIN{f2=f2/100;}NR>1{filter="";if(($5==0 && $6==0)||($8==0 && $9==0)) {filter="COV0";} else {TAF=$6/($5+$6);NAF=$9/($8+$9); if(NAF>f2) filter=filter";MNAF"(100*f2); if(($5+$6)<f4) filter=filter";MRDT"f4; if(($8+$9)<f5) filter=filter";MRDN"f5;}; gsub("^;","",filter); if(filter!="") print $1,$2,$3,$4,filter }' >${PREFIX}.filter1.txt


if [ -s "${PREFIX}.filter1.txt" ]
then
  zcat ${PREFIX}.vcf.gz |sed  '/^##FILTER=<ID=PASS/a ##FILTER=<ID=COV0,Description="No read coverage in the normal sample">\n##FILTER=<ID=MNAF'"$MNAF"',Description="Variant allele frequency in normal higher than cutoff '"$MNAF"'/100">\n##FILTER=<ID=MRDT'"$MRDT"',Description="Read depth in the tumor less than cutoff '"$MRDT"'">\n##FILTER=<ID=MRDN'"$MRDN"',Description="Read depth in the normal less than the cutoff '"$MRDN"'">' |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1$2$3$4]=$5;next;}{if(!/^#/) {key=$1$2$4$5; if(key in a) {gsub("PASS","",$7); if($7=="") $7=a[key]; else $7=$7";"a[key];} else { if($7=="" || $7=="PASS") $7="PASS"; else $7=$7";PASS"} } print $0;}' ${PREFIX}.filter1.txt -  >${PREFIX}.filter1.vcf
else
  zcat ${PREFIX}.vcf.gz |sed  '/^##FILTER=<ID=PASS/a ##FILTER=<ID=COV0,Description="No read coverage in the normal sample">\n##FILTER=<ID=MNAF'"$MNAF"',Description="Variant allele frequency in normal higher than cutoff '"$MNAF"'/100">\n##FILTER=<ID=MRDT'"$MRDT"',Description="Read depth in the tumor less than cutoff '"$MRDT"'">\n##FILTER=<ID=MRDN'"$MRDN"',Description="Read depth in the normal less than the cutoff '"$MRDN"'">' |awk -F "\t" -v OFS="\t" '{if(!/^#/) { if($7=="" || $7=="PASS") $7="PASS"; else $7=$7";PASS"} print $0;}'  >${PREFIX}.filter1.vcf
fi

#zcat ${PREFIX}.vcf.gz |sed  '/^##FILTER=<ID=PASS/a ##FILTER=<ID=COV0,Description="No read coverage in the normal sample">\n##FILTER=<ID=MNAF'"$MNAF"',Description="Variant allele frequency in normal higher than cutoff '"$MNAF"'/100">\n##FILTER=<ID=MRDT'"$MRDT"',Description="Read depth in the tumor less than cutoff '"$MRDT"'">\n##FILTER=<ID=MRDN'"$MRDN"',Description="Read depth in the normal less than the cutoff '"$MRDN"'">' |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1$2$3$4]=$5;next;}{if(!/^#/) {key=$1$2$4$5; if(key in a) {gsub("PASS","",$7); if($7=="") $7=a[key]; else $7=$7";"a[key];} else { if($7=="" || $7=="PASS") $7="PASS"; else $7=$7";PASS"} } print $0;}' ${PREFIX}.filter1.txt -  >${PREFIX}.filter1.vcf 

perl /data/zhangt8/Ref/fpfilter/fpfilter-tool-1.0.1/fpfilter_cram.pl  --vcf-file ${PREFIX}.filter1.vcf  --bam-file ${tumor_bam} --sample $tumor_name --reference $fasta --output ${PREFIX}.filter2.vcf --min-read-pos 0.1  --min-var-freq 0 --min-var-count 4 --min-strandedness 0.01  --max-mm-qualsum-diff 50 --max_var_mm_qualsum 100 --max-mapqual-diff 30 --max-readlen-diff 40 --min-var-dist-3 0.2 

echo "@@ Count the tags: "
bcftools view -H  ${PREFIX}.filter2.vcf |cut -f 7 |sort |uniq -c 

cat ${PREFIX}.filter2.vcf |awk -F "\t" -v OFS="\t" '{gsub("PASS;","",$7); print $0;}' |bcftools view -s ${tumor_name},${normal_name} -O z -o ${PREFIX}.vcf.gz 
bcftools index -f ${PREFIX}.vcf.gz 

rm -rf ${PREFIX}.filter1.vcf ${PREFIX}.filter1.txt ${PREFIX}.filter2.vcf ${PREFIX}.VAF 

#bcftools view -f "PASS" -O z -o ${PREFIX}.filter.vcf.gz  ${PREFIX}.vcf.gz
#bcftools index -f ${PREFIX}.filter.vcf.gz


## annotation steps
convert2annovar.pl -format vcf4old ${PREFIX}.vcf.gz -outfile all.avinput 
zcat ${PREFIX}.vcf.gz | awk -F "\t" -v OFS="\t" 'BEGIN{k=1;}/^#/{print $0; next;}{$3="V"k;k++;print $0;}' >all2.vcf
awk -F "\t" -v OFS="\t"  'BEGIN{k=1;}/^#/{print $0; next;}{$9="V"k;k++;print $0;}' all.avinput >all2.avinput
table_annovar.pl all2.avinput $ANNOVAR_DATA/hg38 --buildver hg38 --outfile all2 --remove --protocol refGene,ensGene,ljb26_all,dbnsfp35c,dbscsnv11,cosmic92_coding,cosmic92_noncoding,exac03,exac03nontcga,1000g2015aug_all,avsnp150,clinvar_20200419,gnomad211_exome,gnomad211_genome,gnomad30_genome  --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --otherinfo  --nastring NA 
head -1  all2.hg38_multianno.txt |transpose -t|nl|awk '$2=="AF" || $2=="ExAC_ALL" || $2=="1000g2015aug_all"{print $1}END{print $1}' |paste -s -d "," |rush ' cut -f {}  all2.hg38_multianno.txt' |awk -F "\t" -v OFS="\t"  'NR>1 && (($1!="NA" && $1>0.001) || ($2!="NA" && $2>0.001) || ($3!="NA" && $3>=0.001) || ($4!="NA" && $4>=0.001)){print $NF}' >MAF_filter.txt 
zcat ${PREFIX}.vcf.gz  >tmp.vcf
bcftools view -H all2.vcf|awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1]=1;next;}($3 in a){ print $0}' MAF_filter.txt - |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1$2$4$5]="MAF001";next;}{if(!/^#/) {key=$1$2$4$5; if(key in a) {if($7=="") $7=a[key];else $7=$7";"a[key]; gsub("PASS;","",$7);}}; print $0;}' - tmp.vcf |sed  '/^##FILTER=<ID=PASS/a ##FILTER=<ID=MAF001,Description="filter by MAF <0.001 in 1000g ExAC gnomAD">' >${PREFIX}.filter3.vcf 
bgzip ${PREFIX}.filter3.vcf 
mv ${PREFIX}.filter3.vcf.gz ${PREFIX}.raw.vcf.gz 
bcftools index -f ${PREFIX}.raw.vcf.gz 
bcftools view -f "PASS" -O z -o ${PREFIX}.vcf.gz ${PREFIX}.raw.vcf.gz
rm -rf ${PREFIX}.vcf.gz.csi 
bcftools index -f ${PREFIX}.vcf.gz

bcftools query -f "%ID\t%FILTER\n" all2.vcf |grep PASS >tmp.pass.txt
awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1]=$0;next;}!($NF in a){print $0}' MAF_filter.txt all2.hg38_multianno.txt | awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next;}($NF in a) || FNR==1 {print $0}' tmp.pass.txt - >${PREFIX}.hg38_multianno_v2.txt 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' all2.vcf |awk -F "\t" -v OFS="\t" '{end=$2+length($3)-1;$2=$2"\t"end; print $0}' |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$NF]=$2;b[$NF]=$3;c[$NF]=$4;d[$NF]=$5;next;}{if(FNR!=1 && (length($4) >1 || $4=="-" || length($5)>1 || $5=="-")) { $2=a[$NF];$3=b[$NF];$4=c[$NF];$5=d[$NF];} print $0}' - ${PREFIX}.hg38_multianno_v2.txt >${PREFIX}.hg38_multianno.txt

zcat ${PREFIX}.raw.vcf.gz | awk -F "\t" -v OFS="\t" 'BEGIN{k=1;}/^#/{print $0; next;}{$3="V"k;k++;if(!/^#/) print $0;}' |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$3]=$7;next;}{if(FNR==1) print $0,"VCFfilter"; else print $0,a[$NF];}' - all2.hg38_multianno.txt >${PREFIX}.hg38_multianno_raw.txt

rm -rf tmp.vcf tmp.pass.txt  all.avinput all2.vcf all2.avinput MAF_filter.txt all2.hg38_multianno.txt 
