#!/bin/bash

module load GATK
module load picard
module load bwa
module load samtools
module load annovar

HAPLOCHECKPATH="~/haplocheckCLI-master"
output_prefix=""
tumor_bam=""
normal_bam=""
tumor_cov=""
original_dir=$( pwd )
if [ -d "/lscratch/$SLURM_JOB_ID" ]
then
	tmp_dir="/lscratch/$SLURM_JOB_ID"
else
	tmp_dir=$original_dir
fi

while getopts "t:n:c:o:" opt; do
	case $opt in
		t) tumor_bam=$OPTARG
		;;
		n) normal_bam=$OPTARG
		;;
		c) tumor_cov=$OPTARG
		;;
		o) out=$OPTARG
		;;
	esac
done

if [[ $tumor_bam == "" ]]
then
	echo "Tumor bam/cram is required"
	exit 1
fi
if [[ $normal_bam == "" ]]
then
	echo "Normal bam/cram is required"
	exit 1
fi

tumor_file="$(awk -F "/" '{print $NF}' <<< $tumor_bam)"
tumor_out_prefix="${tumor_file%.*}"
normal_file="$(awk -F "/" '{print $NF}' <<< $normal_bam)"
normal_out_prefix="${normal_file%.*}"

original_dir=$( pwd )
resource_bundle="/data/Sherlock_Lung/JohnMce/GATK_mito_resources"
hg38="/data/Sherlock_Lung/Share/Reference/hg38/Homo_sapiens_assembly38.fasta"
mito_reference="$resource_bundle/Homo_sapiens_assembly38.chrM.fasta"
shifted_reference="$resource_bundle/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"

realign(){
	input_bam=$1
	output_prefix=$2
	
	#######################################
	## EXTRACT AND PREPARE FOR ALIGNMENT

	# Extract reads aligned to mito chromosome
	if [[ $( awk -F "." '{print $NF}' <<< "$input_bam" ) == "cram" ]];
	then 
	gatk PrintReads \
		-R $hg38 \
		-L chrM \
		--read-filter MateOnSameContigOrNoMappedMateReadFilter \
		-read-filter MateUnmappedAndUnmappedReadFilter \
		-I $input_bam \
		--read-index ${input_bam}.crai \
		-O ${tmp_dir}/${output_prefix}-subsettochrm.bam
	else
	gatk PrintReads \
		-R $hg38 \
		-L chrM \
		--read-filter MateOnSameContigOrNoMappedMateReadFilter \
		-read-filter MateUnmappedAndUnmappedReadFilter \
		-I $input_bam \
		--read-index ${input_bam}.bai \
		-O ${tmp_dir}/${output_prefix}-subsettochrm.bam
	fi

	# Revert sam to unaligned bam
	java -Xmx10g -jar $PICARDJARPATH/picard.jar \
	RevertSam \
	--INPUT ${tmp_dir}/${output_prefix}-subsettochrm.bam \
	--OUTPUT_BY_READGROUP false \
	--OUTPUT ${tmp_dir}/${output_prefix}-revertsam.bam \
	--VALIDATION_STRINGENCY LENIENT \
	--ATTRIBUTE_TO_CLEAR FT \
	--ATTRIBUTE_TO_CLEAR CO \
	--SORT_ORDER queryname \
	--RESTORE_ORIGINAL_QUALITIES false

	# Revert bam to fastq
	java -Xms5000m -jar $PICARDJARPATH/picard.jar \
	SamToFastq \
	--INPUT ${tmp_dir}/${output_prefix}-revertsam.bam \
	--FASTQ ${tmp_dir}/${output_prefix}.fq \
	--INTERLEAVE true \
	--NON_PF true

	rm -r ${tmp_dir}/${output_prefix}_metrics.txt ${tmp_dir}/${output_prefix}-subsettochrm.*
	
	###########################################
	## ALIGN AND CALL TO NORMAL MT GENOME

	# Alignment
	bwa mem -K 100000000 -p -v 3 -t 2 -Y $mito_reference ${tmp_dir}/${output_prefix}.fq | samtools view -bh - > ${tmp_dir}/${output_prefix}-realigned-normal.bam

	bwa_version=$( which bwa | awk -F "/" '{print $(NF-1)}' )

	# Merge info from unaligned bam
	java -Xms3000m -jar $PICARDJARPATH/picard.jar MergeBamAlignment \
	--VALIDATION_STRINGENCY SILENT \
	--EXPECTED_ORIENTATIONS FR \
	--ATTRIBUTES_TO_RETAIN X0 \
	--ATTRIBUTES_TO_REMOVE NM \
	--ATTRIBUTES_TO_REMOVE MD \
	--ALIGNED_BAM ${tmp_dir}/${output_prefix}-realigned-normal.bam \
	--UNMAPPED_BAM ${tmp_dir}/${output_prefix}-revertsam.bam \
	--OUTPUT ${tmp_dir}/${output_prefix}-merged-normal.bam \
	--REFERENCE_SEQUENCE $mito_reference \
	--PAIRED_RUN true \
	--SORT_ORDER "unsorted" \
	--IS_BISULFITE_SEQUENCE false \
	--ALIGNED_READS_ONLY false \
	--CLIP_ADAPTERS false \
	--MAX_RECORDS_IN_RAM 2000000 \
	--ADD_MATE_CIGAR true \
	--MAX_INSERTIONS_OR_DELETIONS -1 \
	--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
	--PROGRAM_RECORD_ID "bwamem" \
	--PROGRAM_GROUP_VERSION "$bwa_version" \
	--PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -v 3 -t 2 -Y $mito_reference ${tmp_dir}/${output_prefix}.fq"  \
	--PROGRAM_GROUP_NAME "bwamem" \
	--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
	--ALIGNER_PROPER_PAIR_FLAGS true \
	--UNMAP_CONTAMINANT_READS true \
	--ADD_PG_TAG_TO_READS false

	rm -r ${tmp_dir}/${output_prefix}-realigned-normal.*

	# Mark Duplicates
	java -Xms4000m -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	--INPUT ${tmp_dir}/${output_prefix}-merged-normal.bam \
	--OUTPUT ${tmp_dir}/${output_prefix}-md-normal.bam \
	--METRICS_FILE ${output_prefix}-normal.metrics \
	--VALIDATION_STRINGENCY SILENT \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
	--ASSUME_SORT_ORDER "queryname" \
	--CLEAR_DT "false" \
	--ADD_PG_TAG_TO_READS false

	rm -r ${tmp_dir}/${output_prefix}-merged-normal.*

	# Sort alignment
	java -Xms4000m -jar $PICARDJARPATH/picard.jar SortSam \
	--INPUT ${tmp_dir}/${output_prefix}-md-normal.bam \
	--OUTPUT ${tmp_dir}/${output_prefix}-normal-sorted.bam \
	--SORT_ORDER "coordinate" \
	--CREATE_INDEX true \
	--MAX_RECORDS_IN_RAM 300000

	rm -r ${tmp_dir}/${output_prefix}-md-normal.*
	cp ${tmp_dir}/${output_prefix}-normal-sorted.* ${original_dir}
	bgzip ${original_dir}/${output_prefix}-normal-sorted.bam

	#############################################
	## ALIGN AND CALL TO SHIFTED MT GENOME

	# Alignment
	bwa mem -K 100000000 -p -v 3 -t 2 -Y $shifted_reference ${tmp_dir}/${output_prefix}.fq | samtools view -bh - > ${tmp_dir}/${output_prefix}-realigned-shifted.bam

	rm -r ${tmp_dir}/${output_prefix}.fq

	bwa_version=$( which bwa | awk -F "/" '{print $(NF-1)}' )

	# Merge info from unaligned bam
	java -Xms3000m -jar $PICARDJARPATH/picard.jar MergeBamAlignment \
	--VALIDATION_STRINGENCY SILENT \
	--EXPECTED_ORIENTATIONS FR \
	--ATTRIBUTES_TO_RETAIN X0 \
	--ATTRIBUTES_TO_REMOVE NM \
	--ATTRIBUTES_TO_REMOVE MD \
	--ALIGNED_BAM ${tmp_dir}/${output_prefix}-realigned-shifted.bam \
	--UNMAPPED_BAM ${tmp_dir}/${output_prefix}-revertsam.bam \
	--OUTPUT ${tmp_dir}/${output_prefix}-merged-shifted.bam \
	--REFERENCE_SEQUENCE $shifted_reference \
	--PAIRED_RUN true \
	--SORT_ORDER "unsorted" \
	--IS_BISULFITE_SEQUENCE false \
	--ALIGNED_READS_ONLY false \
	--CLIP_ADAPTERS false \
	--MAX_RECORDS_IN_RAM 2000000 \
	--ADD_MATE_CIGAR true \
	--MAX_INSERTIONS_OR_DELETIONS -1 \
	--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
	--PROGRAM_RECORD_ID "bwamem" \
	--PROGRAM_GROUP_VERSION "$bwa_version" \
	--PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -v 3 -t 2 -Y $shifted_reference ${tmp_dir}/${output_prefix}.fq"  \
	--PROGRAM_GROUP_NAME "bwamem" \
	--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
	--ALIGNER_PROPER_PAIR_FLAGS true \
	--UNMAP_CONTAMINANT_READS true \
	--ADD_PG_TAG_TO_READS false

	rm -r ${tmp_dir}/${output_prefix}-realigned-shifted.* ${tmp_dir}/${output_prefix}-revertsam.*

	# Mark Duplicates
	java -Xms4000m -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	--INPUT ${tmp_dir}/${output_prefix}-merged-shifted.bam \
	--OUTPUT ${tmp_dir}/${output_prefix}-md-shifted.bam \
	--METRICS_FILE ${output_prefix}-shifted.metrics \
	--VALIDATION_STRINGENCY SILENT \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
	--ASSUME_SORT_ORDER "queryname" \
	--CLEAR_DT "false" \
	--ADD_PG_TAG_TO_READS false

	rm -r ${tmp_dir}/${output_prefix}-merged-shifted.*

	# Sort alignment
	java -Xms4000m -jar $PICARDJARPATH/picard.jar SortSam \
	--INPUT ${tmp_dir}/${output_prefix}-md-shifted.bam \
	--OUTPUT ${tmp_dir}/${output_prefix}-shifted-sorted.bam \
	--SORT_ORDER "coordinate" \
	--CREATE_INDEX true \
	--MAX_RECORDS_IN_RAM 300000

	rm -r ${tmp_dir}/${output_prefix}-md-shifted.*
	cp ${tmp_dir}/${output_prefix}-shifted-sorted.* ${original_dir}
	bgzip ${original_dir}/${output_prefix}-shifted-sorted.bam
}

realign $tumor_bam $tumor_out_prefix
realign $normal_bam $normal_out_prefix
if [[ $out == "" ]]
then
	output_prefix=$tumor_out_prefix
else
	output_prefix=$out
fi

# Calculate median coverage of genome if not specified
if [[ $tumor_cov == "" ]]
then
	read_length=$( samtools view -T $hg38 $tumor_bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | sed -n 500000p | xargs )

	java -Xmx4g -jar $PICARDJARPATH/picard.jar CollectWgsMetrics \
	--INPUT $tumor_bam \
	--VALIDATION_STRINGENCY SILENT \
	--REFERENCE_SEQUENCE $hg38 \
	--READ_LENGTH $read_length \
	--OUTPUT ${tmp_dir}/${output_prefix}_metrics.txt \
	--USE_FAST_ALGORITHM true \
	--INCLUDE_BQ_HISTOGRAM true \
	--THEORETICAL_SENSITIVITY_OUTPUT ${output_prefix}_theoretical_sensitivity.txt

	tumor_cov =$(cat ${tmp_dir}/${output_prefix}_metrics.txt | grep -A 2 "## METRICS" | tail -1 | cut -f 4 )
fi
normal_name=$( samtools view -H ${tmp_dir}/${normal_out_prefix}-normal-sorted.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq )

# Run Mutect 2 outside control region
gatk --java-options "-Xmx10g" Mutect2 \
-R $mito_reference \
-I ${tmp_dir}/${tumor_out_prefix}-normal-sorted.bam \
-I ${tmp_dir}/${normal_out_prefix}-normal-sorted.bam \
--read-filter MateOnSameContigOrNoMappedMateReadFilter \
--read-filter MateUnmappedAndUnmappedReadFilter \
-O ${tmp_dir}/${output_prefix}-normal.vcf \
-L chrM:576-16024 \
--annotation StrandBiasBySample \
--mitochondria-mode \
--normal ${normal_name} \
--max-reads-per-alignment-start 75 \
--max-mnp-distance 0

# Run Mutect 2 within control region (where the circular genome is normally broken to linearize)
gatk --java-options "-Xmx10g" Mutect2 \
-R $shifted_reference \
-I ${tmp_dir}/${tumor_out_prefix}-shifted-sorted.bam \
-I ${tmp_dir}/${normal_out_prefix}-shifted-sorted.bam \
--read-filter MateOnSameContigOrNoMappedMateReadFilter \
--read-filter MateUnmappedAndUnmappedReadFilter \
-O ${tmp_dir}/${output_prefix}-shifted.vcf \
-L chrM:8025-9144 \
--annotation StrandBiasBySample \
--mitochondria-mode \
--normal ${normal_name} \
--max-reads-per-alignment-start 75 \
--max-mnp-distance 0

#################
## MERGE VARIANTS

java -jar $PICARDJARPATH/picard.jar LiftoverVcf \
--I ${tmp_dir}/${output_prefix}-shifted.vcf \
--O ${tmp_dir}/${output_prefix}-shifted-back.vcf \
--R $mito_reference \
--CHAIN $resource_bundle/ShiftBack.chain \
--REJECT ${tmp_dir}/${output_prefix}.rejected.vcf

java -jar $PICARDJARPATH/picard.jar MergeVcfs \
-I ${tmp_dir}/${output_prefix}-shifted-back.vcf \
-I ${tmp_dir}/${output_prefix}-normal.vcf \
-O ${tmp_dir}/${output_prefix}-merged.vcf

gatk MergeMutectStats \
--stats ${tmp_dir}/${output_prefix}-shifted.vcf.stats \
--stats ${tmp_dir}/${output_prefix}-normal.vcf.stats \
-O ${tmp_dir}/${output_prefix}.combined.stats

rm -r ${tmp_dir}/${output_prefix}.rejected.vcf ${tmp_dir}/${output_prefix}-shifted-back.vcf ${tmp_dir}/${output_prefix}-shifted.vcf* ${tmp_dir}/${output_prefix}-normal.vcf*
cp ${tmp_dir}/${output_prefix}-merged.vcf* $original_dir
cp ${tmp_dir}/${output_prefix}.combined.stats $original_dir

###################
## FILTER CALLS

gatk --java-options "-Xmx10g" FilterMutectCalls \
-V ${tmp_dir}/${output_prefix}-merged.vcf \
-R $mito_reference \
-O ${tmp_dir}/${output_prefix}-filtered.vcf \
--stats ${tmp_dir}/${output_prefix}.combined.stats \
--max-alt-allele-count 4 \
--mitochondria-mode \
--min-allele-fraction 0

gatk VariantFiltration \
-V ${tmp_dir}/${output_prefix}-filtered.vcf \
-O ${tmp_dir}/${output_prefix}-first-filter.vcf \
--apply-allele-specific-filters \
--mask $resource_bundle/blacklist_sites.hg38.chrM.bed \
--mask-name "blacklisted_site"

##################################
## CONTAMINATION CHECK
#
mkdir ${tmp_dir}/${output_prefix}-haplocheckCLI

cp ${tmp_dir}/${output_prefix}-first-filter.vcf* ${tmp_dir}/${output_prefix}-haplocheckCLI
cd ${tmp_dir}/${output_prefix}-haplocheckCLI

java -jar $HAPLOCHECKPATH/haplocheckCLI.jar .
sed 's/\"//g' output > output-noquotes
grep "$tumor"  output-noquotes > output-data

if [[ "$( awk '{print $2}' output-data )" == "YES" ]]
then
	contamination_major=$( awk '{print $14}' output-data )
	contamination_minor=$( awk '{print $15}' output-data )
	if (( $(echo "$contamination_major > 0.000" |bc -l) )); then contam=$(echo "1.0 - $contamination_major" | bc -l); else contam=$contamination_minor; fi
	cd $original_dir

	gatk --java-options "-Xmx10g" FilterMutectCalls \
	-V ${tmp_dir}/${output_prefix}-first-filter.vcf \
	-R $mito_reference \
	-O ${tmp_dir}/${output_prefix}-filter-contaminants.vcf \
	--stats ${output_prefix}.combined.stats \
	--max-alt-allele-count 4 \
	--mitochondria-mode \
	--min-allele-fraction 0 \
	--contamination-estimate $contam
fi

cd $original_dir
rm -r ${tmp_dir}/${output_prefix}-haplocheckCLI/*
rm -r ${tmp_dir}/${output_prefix}-haplocheckCLI

######################################
## NORMALIZE AND APPLY LAST FILTERS

if [[ -s ${tmp_dir}/${output_prefix}-filter-contaminants.vcf ]]
then
	gatk LeftAlignAndTrimVariants \
	-R $mito_reference \
	-V ${tmp_dir}/${output_prefix}-filter-contaminants.vcf \
	-O ${tmp_dir}/${output_prefix}-split.vcf \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac
else
	gatk LeftAlignAndTrimVariants \
	-R $mito_reference \
	-V ${tmp_dir}/${output_prefix}-first-filter.vcf \
	-O ${tmp_dir}/${output_prefix}-split.vcf \
	--split-multi-allelics \
	--dont-trim-alleles \
	--keep-original-ac
fi

gatk SelectVariants \
-V ${tmp_dir}/${output_prefix}-split.vcf \
-O ${tmp_dir}/${output_prefix}-splitAndPassOnly.vcf \
--exclude-filtered

cp ${tmp_dir}/${output_prefix}-split.vcf $original_dir
rm -r ${tmp_dir}/${output_prefix}-filter-contaminants.vcf ${tmp_dir}/${output_prefix}-first-filter.vcf

gatk NuMTFilterTool \
-R $mito_reference \
-V ${tmp_dir}/${output_prefix}-splitAndPassOnly.vcf \
-O ${tmp_dir}/${output_prefix}-numt-filtered.vcf \
--autosomal-coverage $tumor_cov

# The original GATK pipeline uses this tool, however it seems to not function properly as of this analysis
# and is thus replaced with the custom python script below
#
#gatk MTLowHeteroplasmyFilterTool \
#-R $mito_reference \
#-V ${tmp_dir}/${output_prefix}-numt-filtered.vcf \
#-O ${output_prefix}-final-call.vcf \
#--max-allowed-low-hets 3

# Argument order is: IN_FILE, VAF_THRESHOLD, MINIMUM_VARS, OUTPUT_PREFIX
python $resource_bundle/MTLowHeteroplasmy.py ${tmp_dir}/${output_prefix}-numt-filtered.vcf 0.01 3 ${tmp_dir}/${output_prefix}

# ANNOVAR steps
#cat ${output_prefix}-final-call.vcf | awk -F "\t" -v OFS="\t" 'BEGIN{k=1;}/^#/{print $0; next;}{$3="V"k;k++;print $0;}' > ${output_prefix}-tmp.vcf
convert2annovar.pl -format vcf4old ${tmp_dir}/${output_prefix}-low-het-filtered.vcf -outfile ${tmp_dir}/${output_prefix}.avinput
awk -F "\t" -v OFS="\t"  'BEGIN{k=1;}/^#/{print $0; next;}{$9="V"k;k++;print $0;}' ${tmp_dir}/${output_prefix}.avinput > ${tmp_dir}/${output_prefix}_2.avinput
annotate_variation.pl --filter --dbtype generic --genericdbfile annovar_gnomad311.txt --buildver hg38 -out ${tmp_dir}/${output_prefix} --otherinfo ${tmp_dir}/${output_prefix}_2.avinput $resource_bundle/

# Argument order is ANNOVAR_IN, VCF_IN, OUTPUT_PREFIX
python $resource_bundle/MTAnnovarFilter.py ${tmp_dir}/${output_prefix}.hg38_generic_dropped ${tmp_dir}/${output_prefix}-low-het-filtered.vcf ${output_prefix}

rm -r ${tmp_dir}/${output_prefix}-numt-filtered.vcf* ${tmp_dir}/${output_prefix}-splitAndPassOnly.vcf* ${tmp_dir}/${output_prefix}-low-het-filtered.vcf
rm -r ${tmp_dir}/${output_prefix}.avinput ${tmp_dir}/${output_prefix}_2.avinput ${tmp_dir}/${output_prefix}.hg38_generic_filtered ${tmp_dir}/${output_prefix}.log

if [[ -d "/lscratch/$SLURM_JOB_ID" ]]
then
	rm -r ${tmp_dir}/${output_prefix}*
fi






