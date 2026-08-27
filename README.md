# we_are_group_12 #SoMosaic

# Exploring False Positives in Somatic Mosaic Structural Variant Calling

This project evaluates and characterizes false positive structural variant (SV)
calls produced by Sniffles2 in long-read sequencing data from the SMaHT Multiple
Individuals in a Mixed Sample (MIMS) benchmark.

The goal is to identify sources of false positive mosaic SV calls and develop
filtering strategies that improve precision while preserving recall.

## Data

The current analysis uses five ONT sequencing datasets aligned to GRCh38 from the MIMS mixture: 
(https://data.smaht.org/data/benchmarking/HapMap)

- 3 datasets from Baylor College of Medicine (BCM) GCC
- 1 dataset from New York Genome Center GCC
- 1 dataset from University of Washington GCC

The datasets were merged before Sniffles2 calling.

The current pilot analysis focuses on chromosome 22 only.

Reference datasets used:

Reference genome: hg38.analysisSet.fa 

Annotate tandem repeats: human_GRCh38_no_alt_analysis_set.trf.bed 
(https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations)


### Benchmark

The MIMS benchmark contains mixtures of six HapMap individuals at defined
abundances to simulate mosaic variation.

We used the **MIMS easy SV benchmark v2** for the primary analysis.
(https://github.com/BCM-HGSC/SMaHT_MIMS/blob/main/README.md)

Files used include:

- `smaht_mims_sv_v2_easy.vcf.gz`
- `smaht_mims_sv_v2_easy.vcf.gz.tbi`
- `smaht_mims_sv_v2_easy_regions.bed`

A GRCh38 no-alt reference genome and its index are also required.


## Workflow
Flowchart draft v3:
<img width="1120" height="475" alt="Screenshot 2026-08-27 at 3 37 43 AM" src="https://github.com/user-attachments/assets/37ff7ba9-e4c2-41c0-8eb4-664d572c0e6d" />

Flowchart draft v2:
<img width="1024" height="628" alt="image" src="https://github.com/user-attachments/assets/672d6a53-a3ff-4ed6-a76b-d59dcf5263d4" />

## Methods
[SoMosaicMethods_1.pdf](https://github.com/user-attachments/files/31518180/SoMosaicMethods_1.pdf)
[SoMosaicMethods_2.pdf](https://github.com/user-attachments/files/31518172/SoMosaicMethods_2.pdf)
[SoMosaicMethods_3.pdf](https://github.com/user-attachments/files/31518171/SoMosaicMethods_3.pdf)



