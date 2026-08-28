# we_are_group_12 #SoMosaic
<img width="1254" height="1254" alt="image" src="https://github.com/user-attachments/assets/441721f8-cb1b-4b1b-97db-4c223b1c9bad" />

# Exploring False Positives in Somatic Mosaic Structural Variant Calling

This project evaluates and characterizes false positive structural variant (SV)
calls produced by Sniffles2 in long-read sequencing data from the SMaHT Multiple
Individuals in a Mixed Sample (MIMS) benchmark.

The goal is to identify sources of false positive mosaic SV calls and develop
filtering strategies that improve precision while preserving recall.

## Quick Start: Core Alignment-Free Streaming Workflow

This repository is organized into an upstream variant extraction workflow and a downstream machine learning classification module located in the sub-directory:

```text
SoMosaic/   (Root)
├── README.md       <-- You are here (Core Overview & Quick Start)
├── LICENSE
├── SoMosaicMethods_20260807.pdf
└── Sniffles2_mosaic_SV_false_positive_filter/  (ML Classification Engine)
    ├── README.md   <-- Dedicated Random Forest guide, feature importance, and scripts
    ├── data
```

### 1. Environment Setup
Ensure your local environment has the required bioinformatics utilities and machine learning packages installed:
```bash
# Install genomic alignment and variant calling tools via conda/mamba
conda install -c biographies bcftools htslib sniffles -y

# Install downstream data matrix and modeling libraries
pip install truvari scikit-learn pandas numpy matplotlib seaborn
```

### 2. SV calling with Sniffles2 Mosaic
Instead of converting CRAM files to intermediate BAM files, Sniffles2 calls somatic variants directly from the sorted CRAM alignments on-the-fly using full multi-threading support:

```bash
sniffles \
  --input /PATH_TO/SMHTHAPMAP6-X-X-NN-D001-uwsc-SMAFIUU2RO1T-sentieon_minimap2_202308.01_GRCh38.aligned.sorted.cram \
  --vcf /PATH_TO/HapMap6.ONT.uwsc-SMAFIUU2RO1T.WGS.sniffles.mosaic.vcf.gz \
  --snf /PATH_TO/HapMap6.ONT.uwsc-SMAFIUU2RO1T.WGS.sniffles.mosaic.snf \
  --sample-id HapMap6_ONT_uwsc_SMAFIUU2RO1T_WGS \
  --tandem-repeats /PATH_TO/human_GRCh38_no_alt_analysis_set.trf.bed \
  --reference /PATH_TO/GRCh38.no_alt_analysis_set.fa \
  --threads 16 \
  --mapq 20 \
  --qc-output-all \
  --output-rnames \
  --mosaic \
  --mosaic-include-germline \
  --mosaic-af-min 0.0
```

### 3. Truth Set Benchmarking (PASS-Only Mode)
Prepare the outputs and evaluate somatic variant calling accuracy metrics against the high-confidence SMaHT MIMS easy v2 intervals utilizing `truvari bench`:

```bash
# Ensure comparison files are block-gzipped and indexed for random access
bcftools sort /PATH_TO/HapMap6.ONT.uwsc-SMAFIUU2RO1T.WGS.sniffles.mosaic.vcf.gz -Oz -o /PATH_TO/HapMap6.ONT.uwsc-SMAFIUU2RO1T.WGS.sniffles.mosaic.bgzf.vcf.gz
tabix -p vcf /PATH_TO/HapMap6.ONT.uwsc-SMAFIUU2RO1T.WGS.sniffles.mosaic.bgzf.vcf.gz

# Execute Truvari evaluation restricted to confidently passing calls
truvari bench \
  -b /PATH_TO/smaht_mims_sv_v2_easy.bgzf.vcf.gz \
  --passonly \
  --includebed /PATH_TO/smaht_mims_sv_v2_easy_regions.bed \
  -c /PATH_TO/HapMap6.ONT.uwsc-SMAFIUU2RO1T.WGS.sniffles.mosaic.bgzf.vcf.gz \
  -o /PATH_TO/truvari_dir_passonly_all
```

### 4. Running the Machine Learning Classifier
Once Truvari populates your True Positive (`tp-call.vcf`) and False Positive (`fp-call.vcf`) outputs in the benchmarking directory, route into our classification subfolder:
```bash
# Move into the subfolder containing feature matrices and models
cd Sniffles2_mosaic_SV_false_positive_filter/

# Run the localized extraction and training scripts detailed in the subfolder README

 **For advanced feature descriptions, training matrices, and Random Forest hyperparameter flags, please proceed to the [Classification Subdirectory README](./Sniffles2_mosaic_SV_false_positive_filter/README.md).**

```

## Pipeline Feature Engineering Matrix

The downstream classification framework constructs tabular dataframes from 23 high-impact features to train the model:

| Feature Dimension | Target Metric Labels | Description / Extraction Logic |
| :--- | :--- | :--- |
| **Call Quality** | `QUAL`, `FMT_GQ` | Overall base quality scores and genotype confidences. |
| **Read Topology** | `SUPPORT`, `SUPPORT_SA`, `SUPPORT_LONG` | Split-alignment counts and soft-clipped long insertion tracking. |
| **Allelic Ratios** | `VAF`, `dr`, `dv`, `depth` | Variant Allele Fraction and directional read distribution mapping. |
| **Alignment Noise**| `NM`, `STDEV_POS`, `STDEV_LEN` | Mean mismatches and boundary start/length coordinate variances. |
| **Coverage Ratios**| `cov_ratio` | Dynamic flanking vector: \(\frac{\text{cov\_center}}{\text{mean}(\text{cov\_up}, \text{cov\_down})}\) |
| **Spatial Density**| `nb_n_1kb` to `nb_n_100kb` | Multi-scale structural crowding counts around the variant point. |
| **Redundancy Tracking**| `nb_best_jaccard`, `nb_best_size_ratio` | Overlap intervals and nested event metrics to catch duplicates. |

---

## Project Performance Summary

Validation tests utilize independent chromosome holdout partitions to ensure model parameters generalize accurately across varying genomic topologies:

| Evaluation Set | Target Cohort | Precision (Default) | Precision (RF Filtered) | Recall (Default) | Recall (RF Filtered) | F1-Score |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Pilot Set** | Chromosome 22 | xx.x% | **xx.x%** | xx.x% | **xx.x%** | **xx.x%** |
| **Holdout Evaluation**| Chromosome 21 | xx.x% | **xx.x%** | xx.x% | **xx.x%** | **xx.x%** |
| **Genome Scale** | Whole Genome (WGS)| xx.x% | **xx.x%** | xx.x% | **xx.x%** | **xx.x%** |

---

## Source Data & Repositories

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

<img width="1116" height="517" alt="Screenshot 2026-08-28 at 1 40 42 PM" src="https://github.com/user-attachments/assets/b5c1a94b-b88e-4a20-bc1f-55635a9d2030" />

## Methods
<img width="6000" height="3375" alt="SoMosaicMethods_1" src="https://github.com/user-attachments/assets/a9a95a8f-ee83-4306-9acc-b619b14a79ef" />
<img width="6000" height="3375" alt="SoMosaicMethods_2" src="https://github.com/user-attachments/assets/8ee26c18-6e8c-4972-82c5-0e577b4e3a12" />
<img width="6000" height="3375" alt="SoMosaicMethods_3" src="https://github.com/user-attachments/assets/0513e10c-60e0-49c7-9de5-7b73aae49c0b" />





