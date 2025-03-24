# Bioinformatics_Projects
# MOFA-based Multi-Omics Integration for Cardiovascular Disease (CVD)

## Overview
This project utilizes the Multi-Omics Factor Analysis (MOFA) framework to integrate and analyze multiomics data related to cardiovascular disease (CVD). MOFA is a powerful tool for identifying hidden factors that drive variability across multiple data modalities, such as genomics, transcriptomics, and proteomics.

## Objectives
- Integrate multiomics datasets (e.g., genomics, transcriptomics, and proteomics) related to CVD.
- Identify key molecular signatures and hidden factors associated with CVD.
- Perform exploratory analysis and visualization of multiomics data.
- Gain insights into disease mechanisms using factor analysis.

## Data Sources
- **Genomics:** Whole genome sequencing (WGS) variant calls - PRJNA264546 with sample id SRR1635389
- **Transcriptomics:** RNA-seq differential expression data (DESeq2 results)- PRJNA394884:  SRR5837852,  SRR5837854,  SRR5837860,  SRR5837866
- **Proteomics:** Mass spectrometry-based quantification (MaxQuant results)- PXD008934

## Workflow
1. **Data Preprocessing**
   - Normalize omics data (e.g., log-transformation, batch correction).
   - Filter features with low variance.
   - Format data into MOFA-compatible structure.
   
2. **MOFA Model Training**
   - Define data views (genomics, transcriptomics, proteomics).
   - Train MOFA model to infer latent factors.
   - Assess model performance and variance explained.
   
3. **Factor Interpretation**
   - Identify key latent factors linked to CVD.
   - Perform gene set enrichment and pathway analysis.
   - Correlate factors with clinical traits.

4. **Visualization & Interpretation**
   - Generate factor heatmaps and scatter plots.
   - Perform hierarchical clustering of latent factors.
   - Create factor association plots with clinical metadata.

## Dependencies
- R (>= 4.0)
- MOFA2 R package
- Bioconductor packages (e.g., limma, DESeq2)
- ggplot2, pheatmap for visualization

## Installation
```r
# Install MOFA2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MOFA2")

# Load MOFA2
library(MOFA2)
```

## Usage
1. **Prepare input data:** Ensure each omics dataset is formatted as a matrix.
2. **Run MOFA model:** Train the model using the provided scripts.
3. **Analyze results:** Interpret the latent factors and perform enrichment analysis.

## Expected Outcomes
- Identification of multiomics-driven molecular signatures in CVD.
- Discovery of novel biomarkers and therapeutic targets.
- Improved understanding of the molecular basis of CVD.


## Contact
For queries, reach out via priyadarshinikilari@gmail.com or GitHub issues.

