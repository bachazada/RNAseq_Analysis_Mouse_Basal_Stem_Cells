# RNAseq Analysis of Basal Stem Cells in Mice

This repository contains code, data, and results for analyzing RNA-seq expression profiles of basal stem-cell enriched cells (B) and committed luminal cells (L) in the mammary gland of virgin, pregnant, and lactating mice. Six groups are present, with one for each cell type and mouse status combination. Each group contains two biological replicates.

## Overview
The data used in this study has already been aligned to the mouse genome, and the featureCounts tool was used to count reads mapped to mouse genes from Refseq annotation.

### **Main Objectives**:
1. Normalize RNA-seq data using edgeR.
2. Visualize expression data and identify composition bias.
3. Perform differential expression analysis between basal pregnant and basal lactating cells.
4. Annotate significant genes with biological information.

## **Project Files**:
### **Scripts**:
- `scripts/analysis.R`: Contains the R script for RNA-seq data analysis, from normalization to differential expression analysis and visualization.
  
### **Data**:
- `data/data_mat.rds`: RNA-seq read counts data for each sample.
- `data/sample_info.rds`: Metadata on the samples, including cell types and experimental conditions.

### **Results**:
- `results/log_cpm_vs_voom.pdf`: Boxplot comparing unnormalized logCPM and voom-transformed logCPM&#8203;:contentReference[oaicite:0]{index=0}.
- `results/M_versus_A_cpm.pdf`: Mean-Difference (M-A) plot for CPM-normalized counts&#8203;:contentReference[oaicite:1]{index=1}.
- `results/M_versus_A_TMM.pdf`: Mean-Difference (M-A) plot after TMM normalization&#8203;:contentReference[oaicite:2]{index=2}.
- `results/limma_result.csv`: Results of the differential expression analysis, including p-values and fold changes.
- `results/significance_plots.pdf`: Volcano and M-A plots of the top 100 most differentially expressed genes&#8203;:contentReference[oaicite:3]{index=3}.

### **Session Information**:
- `session_info.txt`: R session details, including R version and loaded packages for reproducibility&#8203;:contentReference[oaicite:4]{index=4}.

