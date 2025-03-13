## RNA-Seq log2TPM Data (`rna_seq_log2tpm_data.RData`)

This R object contains an RNA-Seq quantification matrix in log2-transformed TPM values (`rna_seq_log2tpm_data`), derived from samples used in this manuscript as part of the **Sherlock-Lung Study**.

### Sample Naming Convention:
The columns (sample names) within the matrix follow the format:

- `RNA-Seq_sample_name`: Unique identifier for the RNA-Seq sample.
- `WGS_sample_name`: Corresponding whole-genome sequencing (WGS) sample identifier.
- `RNA-Seq_sample_type`: Type of the sample (`Tumor`, `Normal`, etc.).

### How to Load and Inspect Data in R:

```r
load("rna_seq_log2tpm_data.RData")

# Inspect the first few rows and columns:
rna_seq_log2tpm_data[1:4, 1:4]

# A tibble: 4 × 4
  Gene     `SC622881:NSLC-0264-T01:Normal` `SC622883:NSLC-0274-T01:Normal` `SC622889:NSLC-0264-T01:Tumor`
  <chr>                              <dbl>                           <dbl>                          <dbl>
1 5S_rRNA                           -3.79                           -3.80                          -3.40 
2 A1BG                              -0.218                          -0.548                         -0.278
3 A1BG-AS1                           1.21                            1.50                           1.81 
4 A2M                               11.4                            11.6                            9.50 
