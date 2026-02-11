# GEO-Differential-Expression

\# Differential Gene Expression Analysis using GEOquery and limma

\## 📌 Introduction

This project demonstrates a bioinformatics workflow using R and Bioconductor packages.  

We analyze RNA-seq data from the Gene Expression Omnibus (GEO) dataset \*\*GSE70970\*\*.  

The goal is to identify differentially expressed genes and visualize them with volcano plots and heatmaps.

---

\## 🛠️ Methods

1\. \*\*Data Retrieval\*\*  

&nbsp;  - Used `GEOquery` to download expression data and sample metadata.  

2\. \*\*Basic Analysis\*\*  

&nbsp;  - Calculated average gene expression.  

&nbsp;  - Visualized distribution with histograms.  

3\. \*\*Differential Expression\*\*  

&nbsp;  - Applied `limma` for linear modeling and empirical Bayes statistics.  

&nbsp;  - Generated a ranked list of differentially expressed genes.  

4\. \*\*Visualization\*\*  

&nbsp;  - \*\*Volcano Plot\*\*: highlights significantly up/downregulated genes.  

&nbsp;  - \*\*Heatmap\*\*: shows clustering of top 50 DE genes across samples.  


---

\## 📊 Results

\- \*\*Volcano Plot\*\*: Genes with |log2FC| > 1 and adjusted p-value < 0.05 are highlighted in red.  

\- \*\*Heatmap\*\*: Clear clustering of samples based on gene expression patterns.  

\### Example Outputs

!\[Volcano Plot](volcano\_plot.png)  

!\[Heatmap](heatmap.png)  

Differential expression results are saved in:

\- `top\_genes.csv`  

---

\## 🚀 How to Run

1\. Clone this repository:

&nbsp;  ```bash

&nbsp;  git clone https://github.com/AmruthaKantheti/GEO-Differential-Expression.git

2\. Open the project in RStudio.

3\. Install required packages:

if (!requireNamespace("BiocManager", quietly = TRUE))

&nbsp;   install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"))

install.packages("pheatmap")

4.Run the script:

source("geo\_analysis.R")

📌 Dependencies
- R (≥ 4.0)
- GEOquery
- limma
- pheatmap

📖 Interpretation
- The volcano plot highlights genes with significant fold changes.
- The heatmap shows clustering of samples, which may correspond to biological conditions (e.g., tumor vs normal).

📂 Project Structure
GEO-Differential-Expression/
│── README.md
│── geo_analysis.R
│── plots/
│   ├── volcano_plot.png
│   └── heatmap.png
│── results/
│   └── top_genes.csv


✨ Future Work
- Perform Gene Ontology (GO) enrichment analysis.
- Explore pathway analysis for significant genes.
- Extend workflow to other GEO datasets.




