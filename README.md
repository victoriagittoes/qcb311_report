# qcb311_report
For this paper (https://doi.org/10.1186/s13059-016-1033-x), characterize the authors’ work, and investigate the data set deeper with your own analysis.

## Stage 1: Defining the DE state
* .R file: final_project.R

### 0. set up data and Seurat object
### 1. QC
- from violin plots, I filtered by: 
### 2. Normalization
- log-normalization (paper did median-by-ratio normalization)
### 3. PCA
- PC4 v. PC6 (change PC axes accordingly)
- SUpplementary Figure: PC1 - PC7 (all combinations)
  - identified that PC4 differentiates DE cells from others
### 4. UMAP
### 5. Hierarchical Clustering
### 6. Supplementary Figures of DE markers found in mine and paper's analysis
- Violin Plots
- UMAPs

## Stage 2: Investigating when DE state emerges
* .R file: pseudotime_monocle.R

### 0. install Monocle if needed
### 1. set up data and Monocle object
### 2. PCA
### 3. Monocle Tutorial
- from https://cole-trapnell-lab.github.io/monocle-release/
### 4. Default Pseudotime
- not included in final report
### 5. Pseudotime for each time point
- time points: 00h, 12h, 24, 36h, 72h, 96h
- not included in final report
### 6. Pseudotime for each time point, and then combine into a singular plot
- genes of interest: POU5F1, T, CXCR4, SOX17 (same as those shown in paper)
### 7. Violin plots of expression of selected genes per time point
