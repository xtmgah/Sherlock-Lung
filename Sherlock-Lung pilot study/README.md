# Sherlock-Lung Pilot Study
Bioinformatic analysis scripts for Sherlock-Lung pilot study. We inlcuded four differnt somatic callers genome wide somatic varaint calling. In summary, the analysis-ready BAM files were processed using four different algorithms, including MuTect, MuTect, Strelka, and TNscope. To improve the performance of the variant calling, we used Sentieon’s genomics package to run MuTect, MuTect and TNscope. Only those SNVs that passed calling by a minimum of three algorithms were kept. To reduce false positive calling, we applied this in-house filtering script similar to our previous publication. To summarize, variant calling was considered only at the genome positions with 1) read depth >12 in tumor and >6 in normal samples; and 2) variant read count >5 in tumor and VAF <0.02 in normal samples. The filtered variants were annotated with Oncotator and ANNOVAR. For the indel calling, only variants called by three algorithms were kept (MuTect2, TNscope, and Strelka). UPS-indel algorithm was used to compare and combine different indel call sets. Similar filtering steps as those used for SNV calling were also applied to indel call. We recommend to remove possible germline variants from the called somatic variants, using the dbSNP138, 1000 genomes, ExAC, gnomAD database, or any in-house germline variant database for commonly occurring SNPs (somatic variant frequency<0.001). 


## Examples

### for SNVs
sh variant-filter-snv.sh /PATH_to_TNSNV_VCF/ /PATH_to_TNHAPLOTYPER_VCF/ /PATH_to_TNSCOPE_VCF/ /PATH_to_STRELKA_SNVs_VCF/ /PATH_to_STRELKA_INDEL_VCF/ /PATH_to_NORMAL_VCF/ /PATH_to_Tumor_BAM/ OUTPUT_PREFEIX

sh variant-filter-snv.sh /data/zhangt8/NSLC/Calling/Sentieon/TNseq/NSLC-0002-T01/NSLC-0002-T01-tnsnv.vcf.gz /data/zhangt8/NSLC/Calling/Sentieon/TNseq/NSLC-0002-T01/NSLC-0002-T01-tnhaplotyper.vcf.gz /data/zhangt8/NSLC/Calling/Sentieon/TNseq/NSLC-0002-T01/NSLC-0002-T01-tnscope.vcf.gz /data/zhangt8/NSLC/Calling/strelka/NSLC-0002-T01/strelka_output/results/variants/somatic.snvs.vcf.gz /data/zhangt8/NSLC/Calling/strelka/NSLC-0002-T01/strelka_output/results/variants/somatic.indels.vcf.gz /data/ITEB_Lung_WGS/Rawdata/Wave1/SC349575.bam /data/ITEB_Lung_WGS/Rawdata/Wave1/SC349706.bam NSLC-0002-T01

## for INDELs
sh variant-filter-snv.sh /PATH_to_TNHAPLOTYPER_VCF/ /PATH_to_TNSCOPE_VCF/ /PATH_to_STRELKA_INDEL_VCF/ /PATH_to_Tumor_BAM/ /PATH_to_NORMAL_VCF/ OUTPUT_PREFEIX

sh variant-filter-indel.sh /data/zhangt8/NSLC/Calling/Sentieon/TNseq/NSLC-0123-T01/NSLC-0123-T01-tnhaplotyper.vcf.gz /data/zhangt8/NSLC/Calling/Sentieon/TNseq/NSLC-0123-T01/NSLC-0123-T01-tnscope.vcf.gz /data/zhangt8/NSLC/Calling/strelka/NSLC-0123-T01/strelka_output/results/variants/somatic.indels.vcf.gz /data/ITEB_Lung_WGS/Rawdata/Wave1/SC349808.bam /data/ITEB_Lung_WGS/Rawdata/Wave1/SC349821.bam NSLC-0123-T01
