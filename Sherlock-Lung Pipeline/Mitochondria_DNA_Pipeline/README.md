# Mitochondria Pipeline

Adapted from the GATK mitocondria variant calling pipeline, found here: https://app.terra.bio/#workspaces/help-gatk/Mitochondria-SNPs-Indels-hg38, for use on Biowulf. Additionally, further scripts were added to filter out known SNPs present in ANNOVAR and annotate variants with Funcotator.

## Usage
Dependencies: GATK 4.2.0 or higher, haplocheckCLI, ANNOVAR 

Arguments:
-t: path to the tumor alignment file, required

-n: path to the normal alignment file, required

-c: autosomal coverage for the tumor alignment. If not specified this will be calculated using Picard CollectWgsMetrics but note that this takes quite a long time and will account for much of the pipeline's runtime.

-o: output prefix for all outputs from the pipeline
