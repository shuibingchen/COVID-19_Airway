# COVID-19_Airway

This repository contains the R scripts necessary to perform the analysis in our
paper [Duan, Xiaohua et al. "An Airway Organoid-Based Screen Identifies a Role
for the HIF1a-Glycolysis Axis in SARS-CoV-2
Infection." *Cell Reports*, in press](https://github.com/shuibingchen/COVID-19_Airway), as described in the supplementary methods and main text.

### Input data
The single cell RNA-seq data were generated with the 10X Chromium and
pre-processed using 10X cellranger pipeline. The raw data are available in the
GEO database with
accession#[GSE160231](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160231).

### Requirements
The following R packages are needed:
- Seurat
- scran
- scater
- future
- dplyr
- reshape2
- magrittr
- ggplot2
- pheatmap
- cowplot
- RColorBrewer

