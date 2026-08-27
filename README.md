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
<img width="6000" height="3375" alt="SoMosaicMethods_3" src="https://github.com/user-attachments/assets/0513e10c-60e0-49c7-9de5-7b73aae49c0b" />
<img width="6000" height="3375" alt="SoMosaicMethods_2" src="https://github.com/user-attachments/assets/8ee26c18-6e8c-4972-82c5-0e577b4e3a12" />
<img width="6000" height="3375" alt="SoMosaicMethods_1" src="https://github.com/user-attachments/assets/a9a95a8f-ee83-4306-9acc-b619b14a79ef" />





