# WGS mutation calling pipeline for Sherlock-Lung Study (N=1217; Genome build: GRCh38)
In summary, we employed four different somatic callers for genome-wide somatic variant calling. Specifically, the analysis-ready BAM files were processed using four different algorithms: MuTect, MuTect, Strelka, and TNscope. To improve the performance of the variant calling, we used Sentieon’s genomics package to run MuTect, MuTect, and TNscope. Only SNVs that were called by at least three algorithms were retained to reduce false positives. To further reduce false positives, we applied an in-house filtering script similar to our previous publications. Overall, variant calling was performed only at genome positions that met the following criteria: 1) read depth greater than 12 in tumor and greater than 6 in normal samples; and 2) variant read count greater than 5 in tumor and variant allele frequency less than 0.02 in normal samples. The filtered variants were annotated with ANNOVAR. For indel calling, we retained only those variants that were called by three algorithms (MuTect2, TNscope, and Strelka). We used the UPS-indel algorithm to compare and merge different indel call sets. We applied similar filtering steps to indel calling as those used for SNV calling. We recommend removing potential germline variants from the called somatic variants using the dbSNP, 1000 Genomes, ExAC, gnomAD database, or any in-house germline variant database for commonly occurring SNPs (somatic variant frequency < 0.001).

## Examples

### for both SNVs and INDELs

sh variant-filter-snv.sh /PATH_to_TNSNV_VCF/ /PATH_to_TNHAPLOTYPER_VCF/ /PATH_to_TNSCOPE_VCF/ /PATH_to_STRELKA_SNVs_VCF/  /PATH_to_STRELKA_INDELs_VCF/ /PATH_to_NORMAL_CRAM/ /PATH_to_Tumor_CRAM/ OUTPUT_PREFEIX


For example, 

sh variant-filter-multiple.sh 
/data/zhangt8/NSLC2/Calling/Sentieon/TNseq/NSLC-0269-T01/NSLC-0269-T01-tnsnv.vcf.gz 
/data/zhangt8/NSLC2/Calling/Sentieon/TNseq/NSLC-0269-T01/NSLC-0269-T01-tnhaplotyper.vcf.gz
/data/zhangt8/NSLC2/Calling/Sentieon/TNseq/NSLC-0269-T01/NSLC-0269-T01-tnscope.vcf.gz
/data/zhangt8/NSLC2/Calling/strelka/NSLC-0269-T01/strelka_output/results/variants/somatic.snvs.vcf.gz
/data/zhangt8/NSLC2/Calling/strelka/NSLC-0269-T01/strelka_output/results/variants/somatic.indels.vcf.gz
/data/ITEB_Lung_WGS/FireCloud_downloads/Paper2/NCIContract_DCEG_NSLC_84samples_May2019/NSLC-AI89-NT1-A-1-1-D-A77N-36/v1/NSLC-AI89-NT1-A-1-1-D-A77N-36.cram /data/ITEB_Lung_WGS/FireCloud_downloads/Paper2/NCIContract_DCEG_NSLC_84samples_May2019/NSLC-AI89-TTP1-A-1-1-D-A77N-36/v1/NSLC-AI89-TTP1-A-1-1-D-A77N-36.cram
NSLC-0269-T01



