
## Deciphering lung adenocarcinoma evolution and the role of LINE-1 retrotransposition

This repository contains R code used to generate the visualizations for all main figures developed as part of the Sherlock-Lung project. The corresponding R scripts also include additional plots for relevant supplementary figures.

For reproducibility purposes, we have also included the R Markdown (.Rmd) files along with their rendered HTML outputs for each script.

**To see the full display and formatting correctly, please download the HTML file and open it on your local computer.**

## Figures included ## 

**Fig. 1: Evolutionary dynamics of lung cancer (Fig1_Tumor_Evolution.R)** 
a) Sankey diagram illustrating high-clonality WGS data summary from Sherlock-Lung.
b) Proportion of tumor samples exhibiting whole genome doubling (WGD) across AS_N, EU_N, EU_S and “Others”.
c) Distribution of the percentage of mutations with a copy number 2, considering only mutations attributed to clock-like signatures (SBS1 and SBS5).
d) Box plots depicting the proportion of total mutations attributed to different clonal statuses.
e) Enrichment of early clonal mutations within driver genes in never-smokers and smokers.
f) Evolutionary models displaying the recurrent temporal order of driver genes, as inferred by the ASCETIC framework.
g) Dynamic mutational processes during clonal and subclonal tumor evolution.

**Fig. 2: Features associated with lung tumor latency (Fig2_Tumor_Latency_related_analyses.R)**
a) Associations between tumor latency and driver gene mutation status.
b) Forest plot illustrating associations between tumor latency and EGFR mutation status, adjusted for sex and tumor purity (percentage of cancer cells within a tumor sample).
c) Box plots displaying estimated tumor latency separated by EGFR mutation status and sex.
d) Associations between tumor latency and the presence of specific mutational signatures.
e) A multivariate regression analysis examining the relationship between tumor latency and various factors, including sex, ancestry, smoking, EGFR mutations, KRAS mutations, and mutational signature ID2.

**Fig. 3: Characterization of tumors with mutational signature ID2 (Fig3_ID2_characterization.R)**
a) Relationships between mutational signature ID2 presence and gene expression of tumor proliferation markers, analyzed from tumor and normal tissue RNA-Seq data.
b) Pearson correlation between the number of deletions attributed to mutational signature ID2 and the gene expression of tumor proliferation markers.
c) Kaplan-Meier survival curves for overall survival, stratified by the presence of mutational signature ID2.
d) Enrichment of tumor metastasis in tumors with mutational signature ID2.
e) Enrichment of genomic alterations (WGD=whole genome doubling; SCNA=somatic copy number alterations) in tumors with mutational signature ID2, determined through logistic regression and adjusted for ancestry, sex, smoking status, age, and tumor purity.
f) Distribution of estimated hypoxia scores between tumors with ID2 signatures and those without.

**Fig. 4: Association between L1 retrotransposition and mutational signatures ID2 and ID1 (Fig4_LINE1_detection_association.R)**
a) This panel illustrates a sample (NSLC-0622-T01) as an example of a tumor harboring L1 insertions from a germline source.
b) This section presents another example (tumor sample NSLC-0832-T01) predominantly characterized by somatic-source L1 insertions.
c) Distribution of retrotransposable sources of L1 insertions.
d) Enrichment of mutational signatures in tumors with germline source L1 insertions.
e) Pearson correlation between deletions attributed to signature ID2 and insertions attributed to signature ID1.
f-g) Mutational signature profiles and motifs for mutational signature ID1 and ID2, respectively.
h) Pearson correlation between deletions attributed to ID2 and total somatic L1 insertions.

**Fig. 5: Activation of germline L1 retrotransposition due to DNA demethylation of L1 promoter (Fig5_LINE1_Methylation.R)**
a) Diagram illustrating the transposon mobilization mechanism for long interspersed element 1 (L1).
b) Validation of DNA methylation levels in the promoter region of germline L1 insertions from germline source in chr22q12.1, conducted via targeted bisulfite sequencing.
c) Box plot shows DNA median methylation levels across the genome locations of the CpG island on chr22q12.1, stratified by the sample type and ID2 status.
d) Total L1 RNA expression estimated from RNA-Seq data, stratified by sample type and group.
e) Total L1 RNA expression differs between tumors with and without the ID2 signature.

**Fig. 6: ZNF695 upregulation in tumors and its association with mutational signature ID2 (Fig6_KZFPs_related_analyses.R)**
a) Analysis of differentially expressed KZFP protein coding genes between tumors with ID2 signatures and those without.
b) Box plots illustrate the differential expression of ZNF695 among normal tissue or blood, tumors without ID2 signatures, and tumors with ID2 signatures.
c) Pearson correlations between KZFP protein-coding gene expression and indels attributed to mutational signature ID2.
d) Correlation between ZNF695 expression and deletions attributed to mutational signature ID2.
e) Correlation between ZNF695 RNA-Seq expression and the median of DNA methylation levels across the genome locations of the CpG island on chr22q12.1.
f) Differentially expressed ZNF695 target genes identified between tumors with and without mutational signature ID2.

---

## Instructions for Reproducibility ##

All R scripts are written in a generic way to allow users to reproduce the figures presented in this tumor evolution paper. To run the code successfully, please follow the folder structure and guidelines below:

Input Data: Place all required input data files in the data/ directory.

Output Files: All visualizations and generated plots will be saved in the output/ directory.

Functions: The functions/ folder contains additional helper functions required for visualization. These will be sourced automatically by the scripts.

Note:
The data/ folder is expected to contain most of the whole-genome sequencing (WGS) analysis results from the Sherlock-Lung study. These are large files and are not included in this repository due to data access restrictions. To obtain the raw WGS data, users must apply through dbGaP (application details are provided in our manuscript). All sequencing raw data have been deposited in dbGaP in accordance with NIH guidelines and data use restrictions. For any questions, please contact the Principal Investigator (Dr. Maria Teresa Landi) of the Sherlock-Lung study directly.

As an example, we have included the normalized RNA-Seq data in R object format, as shown below.

#### RNA-Seq log2TPM Data (`rna_seq_log2cpm_data.RData`)

This R object contains an RNA-Seq quantification matrix in log2-transformed CPM values (`rna_seq_log2cpm_data`), derived from samples used in this manuscript as part of the **Sherlock-Lung Study**.

#### Sample Naming Convention:
The columns (sample names) within the matrix follow the format:

- `RNA-Seq_sample_name`: Unique identifier for the RNA-Seq sample.
- `WGS_sample_name`: Corresponding whole-genome sequencing (WGS) sample identifier.
- `RNA-Seq_sample_type`: Type of the sample (`Tumor`, `Normal`, etc.).

#### How to Load and Inspect Data in R:

```r
load("rna_seq_log2cpm_data.RData")

# Inspect the first few rows and columns:
rna_seq_log2cpm_data[1:4, 1:4]

# A tibble: 4 × 4
  Gene     `SC622881:NSLC-0264-T01:Normal` `SC622883:NSLC-0274-T01:Normal` `SC622889:NSLC-0264-T01:Tumor`
  <chr>                              <dbl>                           <dbl>                          <dbl>
1 5S_rRNA                           -3.79                           -3.80                          -3.40 
2 A1BG                              -0.218                          -0.548                         -0.278
3 A1BG-AS1                           1.21                            1.50                           1.81 
4 A2M                               11.4                            11.6                            9.50 
