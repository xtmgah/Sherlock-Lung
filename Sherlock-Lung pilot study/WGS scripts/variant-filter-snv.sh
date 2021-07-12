#!/bin/bash
set -e 
##### Software Variant #######
#module load samtools
#module load GATK/3.8-0 
#module load oncotator/1.9.1.0
#module load snpEff 
##### Database Variant #####
Transcript_override_lists=/data/zhangt8/Ref/Oncotator/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt
##### Input ########

mutect_VCF=$1
mutect2_VCF=$2
tnscope_VCF=$3
strelka_snv_VCF=$4
strelka_indel_VCF=$5
normal_bam=$6
tumor_bam=$7
PREFIX=$8

fasta="/data/zhangt8/APOBEC/Ref/Homo_sapiens_assembly19.fasta"

DIR=$(pwd)/Result/${PREFIX}
MAFOUTPUT=${DIR}/$PREFIX}.maf
SGE_LOG=${DIR}/${PREFIX}.log

mkdir -p $DIR

cd $DIR

### index vcf files #### 
#bcftools index -t $mutect_VCF
#bcftools index -t $mutect2_VCF
#bcftools index -t $tnscope_VCF



#### extract normal and tumor name ####

normal_name=`samtools view -H ${normal_bam} |grep "@RG" |sed 's/.*SM://'|awk '{print $1}' |uniq`
tumor_name=`samtools view -H  ${tumor_bam} |grep "@RG" |sed 's/.*SM://'|awk '{print $1}' |uniq `

### process snp call ######

echo "NORMAL $normal_name" >newname.txt
echo "TUMOR $tumor_name" >>newname.txt 

bcftools view -i 'SVTYPE="."' $tnscope_VCF -O z -o tnscope_VCF.vcf.gz 
bcftools index -t tnscope_VCF.vcf.gz 

bcftools reheader -s newname.txt $strelka_snv_VCF |bcftools view -O z >strelka.somatic.snvs.vcf.gz
bcftools index -t strelka.somatic.snvs.vcf.gz 
java -Xmx4g -jar  $GATK_JAR  -T CombineVariants -R $fasta -genotypeMergeOptions PRIORITIZE --variant:mutect2 $mutect2_VCF --variant:mutect $mutect_VCF  --variant:tnscope tnscope_VCF.vcf.gz --variant:strelka strelka.somatic.snvs.vcf.gz -o combine_snv.vcf -priority mutect2,mutect,tnscope,strelka --minimumN 3 --filteredrecordsmergetype KEEP_IF_ALL_UNFILTERED
bcftools view  -m2 -M2 -i ' set="Intersection" | set="mutect2-mutect-tnscope" |set="mutect2-mutect-strelka" | set="mutect2-tnscope-strelka" |set="mutect-tnscope-strelka"  |set="filterInmutect2-mutect-tnscope-strelka" |set="mutect2-filterInmutect-tnscope-strelka"|set="mutect2-mutect-filterIntnscope-strelka"|set="mutect2-mutect-tnscope-filterInstrelka" ' combine_snv.vcf  >combine_snv_filter.vcf 

#bcftools reheader -s newname.txt $strelka_indel_VCF |bcftools view -O z >strelka.somatic.indels.vcf.gz
#bcftools index -t strelka.somatic.indels.vcf.gz 
#java -Xmx4g -jar  $GATK_JAR  -T CombineVariants -R $fasta -genotypeMergeOptions PRIORITIZE --variant:mutect2 $mutect2_VCF  --variant:tnscope tnscope_VCF.vcf.gz --variant:strelka strelka.somatic.indels.vcf.gz -o combine_indel.vcf -priority mutect2,tnscope,strelka --minimumN 3 --filteredrecordsmergetype KEEP_IF_ALL_UNFILTERED
#bcftools view  -m2 -M2 -i ' set="Intersection" ' combine_indel.vcf >combine_indel_filter.vcf


python /data/zhangt8/Squamous_Cell_Carcinoma/SNV/add_to_vcf_header.py combine_snv_filter.vcf '##INFO=<ID=OVERLAP,Number=0,Type=Flag,Description="Somatic indel possibly overlaps a second indel.">' >combine_snv_filter.vcf2
mv combine_snv_filter.vcf2 combine_snv_filter.vcf

#java -cp $GATK_JAR org.broadinstitute.gatk.tools.CatVariants  -R $fasta -V combine_snv_filter.vcf  -V combine_indel_filter.vcf -assumeSorted -out combine_all_filter.vcf 
mv combine_snv_filter.vcf combine_all_filter.vcf  

bcftools sort combine_all_filter.vcf  -o ${PREFIX}.vcf.gz  -O z 
bcftools index ${PREFIX}.vcf.gz

rm -rf strelka.somatic.snvs.vcf.gz* combine_snv.vcf* combine_snv_filter.vcf* strelka.somatic.indels.vcf.gz* combine_indel.vcf* combine_indel_filter.vcf* combine_all_filter.vcf* tnscope_VCF.vcf.gz*

### filtered by PASS ###

fname1=`bcftools view -h ${PREFIX}.vcf.gz |tail -1|cut -f 10` 
fname2=`bcftools view -h ${PREFIX}.vcf.gz |tail -1|cut -f 11`

if [ "$fname1" != "$tumor_name" ]
then
	echo "Tumor ID in vcf file is not the first"
	java -jar $SNPSIFT_JAR extractFields  ${PREFIX}.vcf.gz CHROM POS REF ALT GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP GEN[0].AD[0] GEN[0].AD[1] GEN[0].DP >${PREFIX}.VAF
	#exit
else 
	echo "Tumor ID in vcf file is the first"
	java -jar $SNPSIFT_JAR extractFields  ${PREFIX}.vcf.gz CHROM POS REF ALT GEN[0].AD[0] GEN[0].AD[1] GEN[0].DP  GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP >${PREFIX}.VAF
fi 

cat ${PREFIX}.VAF|awk -F "\t" -v OFS="\t" 'NR>1{if(($5==0 && $6==0)||($8==0 && $9==0)) print $1,$2,$3,$4,"0","0"; else {TAF=$6/($5+$6);NAF=$9/($8+$9); if(NAF>0.02||$6<5||($5+$6)<12||($8+$9)<6) print $1,$2,$3,$4,TAF,NAF}}' >${PREFIX}.VAF.bad

## TAF<0.07| 

### #### filtered by germline ##### 

bcftools query -f '%CHROM\t%POS\n'  ${PREFIX}.vcf.gz |sort -k1,1 -k2,2n >${PREFIX}.regions.txt
bcftools query -f 'chr%CHROM\t%POS\n'  ${PREFIX}.vcf.gz |sort -k1,1 -k2,2n >${PREFIX}.regions_chr.txt

bcftools query -R ${PREFIX}.regions.txt -f  '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' /data/zhangt8/EAGLE-WGS-Kidney/Germline/full_callset.ICGC.vcf.bgz |awk -F "\t" -v OFS="\t" '{print $0,"ICGC"}' >${PREFIX}.germline_mutaitons.txt
bcftools query -R ${PREFIX}.regions.txt -f  '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' /data/zhangt8/EAGLE-WGS-Kidney/Germline/full_callset.TCGA.vcf.bgz |awk -F "\t" -v OFS="\t" '{print $0,"TCGA"}' >>${PREFIX}.germline_mutaitons.txt
bcftools query -R ${PREFIX}.regions_chr.txt -f  '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' /data/zhangt8/EAGLE-WGS-Kidney/Germline/variants_annotated.vcf.gz |awk -F "\t" -v OFS="\t" '{print $0,"EAGLE"}'|sed 's/^chr//g' >>${PREFIX}.germline_mutaitons.txt
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${PREFIX}.vcf.gz |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1":"$2]=$3"\t"$4;next}{print $0,a[$1":"$2]}' - ${PREFIX}.germline_mutaitons.txt >${PREFIX}.germline_mutaitons.txt2

##### processing snps ######

bcftools view -s $tumor_name ${PREFIX}.vcf.gz >${PREFIX}.vcf

perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${PREFIX}.vcf  >${PREFIX}.var
bam-readcount -q1 -b15 -w1 -l ${PREFIX}.var -f $fasta ${tumor_bam}  >${PREFIX}.readcount
fpfilter.pl --var-file ${PREFIX}.vcf --readcount-file  ${PREFIX}.readcount --output-file ${PREFIX}.fpfilter

awk -F "\t" -v OFS="\t" 'NR==FNR{if($8=="PASS") a[$1":"$2":"$3":"$4]=1;next}{if(/^#/) print $0;else { if(($1":"$2":"$4":"$5 in a)|| length($4)!=1 || length($5)!=1) print $0 }}' ${PREFIX}.fpfilter  ${PREFIX}.vcf | awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1":"$2":"$3":"$4]=1;next}{if(!($1":"$2":"$4":"$5 in a)) print $0 }' ${PREFIX}.VAF.bad  - |awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1":"$2":"$3]=1;next;}{if(/^#/) print $0; else { if(!($1":"$2":"$4 in a)) print $0 }}' ${PREFIX}.germline_mutaitons.txt - >${PREFIX}.filter.vcf  
#awk -F "\t" -v OFS="\t" 'NR==FNR{if($8=="PASS") a[$1":"$2":"$3":"$4]=1;next}{if(/^#/) print $0;else { if(($1":"$2":"$4":"$5 in a)|| length($4)!=1 || length($5)!=1) print $0 }}' ${PREFIX}.fpfilter  ${PREFIX}.vcf | awk -F "\t" -v OFS="\t" 'NR==FNR{a[$1":"$2":"$3":"$4]=1;next}{if(!($1":"$2":"$4":"$5 in a)) print $0 }' ${PREFIX}.VAF.bad  -  >${PREFIX}.filter.vcf 

#oncotator -v --db-dir=$ONCOTATOR_DATASOURCE --input_format=VCF  --canonical-tx-file $Transcript_override_lists --log_name ${PREFIX}.log --override_config ../override_config --annotate-manual="Tumor_Sample_Barcode:${PREFIX}" --annotate-manual="Matched_Norm_Sample_Barcode:${NORMAL}"  ${PREFIX}.filter.vcf ${PREFIX}.maf hg19
#oncotator -v --db-dir=$ONCOTATOR_DATASOURCE --input_format=VCF --tx-mode EFFECT --log_name ${PREFIX}.log --override_config ../override_config --annotate-manual="Tumor_Sample_Barcode:${PREFIX}" --annotate-manual="Matched_Norm_Sample_Barcode:${NORMAL}"  ${PREFIX}.filter.vcf ${PREFIX}.maf2 hg19

rm ${PREFIX}.readcount 

convert2annovar.pl -format vcf4 ${PREFIX}.filter.vcf  -outfile ${PREFIX}.avinput

table_annovar.pl ${PREFIX}.avinput $ANNOVAR_DATA/hg19 \
  --tempdir /lscratch/$SLURM_JOB_ID \
  --thread 8 \
  --buildver hg19 \
  --outfile ${PREFIX} \
  --remove \
  --protocol refGene,ensGene,avsift,ljb26_all,dbnsfp30a,cg46,dbscsnv11,cosmic64,cosmic70,exac03,exac03nontcga,1000g2015aug_all,1000g2012apr_all,snp138,avsnp147,clinvar_20160302,gnomad_exome,gnomad_genome \
  --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
  --nastring ''
