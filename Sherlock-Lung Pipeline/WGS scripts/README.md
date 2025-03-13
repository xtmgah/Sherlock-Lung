# Normal-Tumor Matched WGS Mutation Calling Pipeline
### Sherlock-Lung Study (Genome Build: GRCh38)
This repository provides scripts for preprocessing, mutation calling, and annotation of Whole Genome Sequencing (WGS) data for tumor-normal matched samples.

## Overview
We used four distinct somatic variant callers to identify genome-wide somatic mutations. The analysis-ready BAM files were processed using four algorithms: **MuTect**, **MuTect2**, **Strelka**, and **TNscope**. To enhance calling accuracy, we combined and filtered these variant calls to minimize false positives.

The pipeline consists of three main steps:

- **Step 1 (Preprocessing)**: Converts BAM/CRAM (or raw FASTQ) into analysis-ready BAM.
- **Step 2: Mutation Calling** for both somatic and germline mutations using Sentieon.
- **Step 3: Merging and annotation** of variant calls using multiple callers (Sentieon and Strelka).

---

## Requirements

- **Sentieon software** (license required)
- **Samtools**, **bcftools**, **GATK**, **Strelka**
- **Python** and custom scripts for annotations
- **Human reference genome**: hg38 (GRCh38)

---

## Workflow Overview

### Step 1: Preprocessing Raw Data, Generating Analysis-ready BAM

`pipeline-hg38CRAM-ztw.sh` preprocesses public or raw sequencing data (BAM/CRAM or FASTQ) into analysis-ready BAM files.

### Usage:

```bash
bash pipeline-hg38CRAM-ztw.sh <bam/cram_file> <output_folder> <nthreads> [shortReadremove]

```
<input_bam_or_cram>: Input BAM or CRAM file.

<output_folder>: Directory to save outputs.

'<'nthreads'>': Number of CPU threads.

[shortReadremove]: Optional; set to true to remove reads shorter than 30 bp.

### Example:
```bash
bash pipeline-hg38CRAM-ztw.sh sample.cram SampleOutput 16 true
```

---
### Step 2: Somatic and Germline Mutation Calling
Perform mutation calling using Sentieon algorithms (**MuTect**, **MuTect2**, **TNscope**, etc.):

```bash
bash pipeline-ztw-tumor_normal.sh \
  <tumor_sample_name> <tumor_bam> \
  <normal_sample_name> <normal_bam> \
  <output_prefix> [optional_rgfile]
```

<tumor_sample_name>: Tumor sample ID.

<tumor_bam>: Tumor analysis-ready BAM file from step 1.

<normal_sample_name>: Normal sample ID.

<normal_bam>: Normal BAM file.

<output_prefix>: Prefix to organize output results.

[optional_rgfile]: Optional read group information file if required.

### Example: 
```bash
bash pipeline-ztw-tumor_normal.sh \
  NSLC-0269-T01 tumor.bam \
  NSLC-0269-N01 normal.bam \
  NSLC-0269
```
Outputs include somatic SNVs, indels, structural, and germline variants.

---
### Step 3: Merge and Annotate Mutation Calls
Merge, filter, and annotate variants from multiple callers:

```bash
bash variant-filter-multiple.sh \
  <TNSNV_VCF> \
  <TNHAPLOTYPER_VCF> \
  <TNSCOPE_VCF> \
  <STRELKA_SNVs_VCF> \
  <STRELKA_INDELs_VCF> \
  <NORMAL_CRAM> \
  <TUMOR_CRAM> \
  <OUTPUT_PREFIX>
```
### Examples:
```bash
bash variant-filter-multiple.sh \
  /data/NSLC2/Calling/Sentieon/TNseq/NSLC-0269-T01/NSLC-0269-T01-tnsnv.vcf.gz \
  /data/NSLC2/Calling/Sentieon/TNseq/NSLC-0269-T01/NSLC-0269-T01-tnhaplotyper.vcf.gz \
  /data/NSLC2/Calling/Sentieon/TNseq/NSLC-0269-T01/NSLC-0269-T01-tnscope.vcf.gz \
  /data/NSLC2/Calling/strelka/NSLC-0269-T01/strelka_output/results/variants/somatic.snvs.vcf.gz \
  /data/NSLC2/Calling/strelka/NSLC-0269-T01/strelka_output/results/variants/somatic.indels.vcf.gz \
  /data/NCIContract_NSLC_84samples/NSLC-AI89-NT1-A-1-1-D-A77N-36.cram \
  /data/NCIContract_NSLC_84samples/NSLC-AI89-TTP1-A-1-1-D-A77N-36.cram \
  NSLC-0269-T01
```
Filtered variants will be annotated using ANNOVAR. Final variants include annotations of clinical relevance, frequency, and functional impact.


