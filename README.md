# -bioinfo-internship-log
My internship scripts and notes in bioinformatics lab

### 7.8 - Density Plot & Outlier Detection（DistributionΔ(a-b)|piechart_celltype）
- Read large .txt file with 79,000+ rows
- Calculated difference between two columns
- Removed outliers using quantile and 3-sigma rules
- Plotted density graph using ggplot2
- Learned data structure inspection with summary(), str(), etc.

### 7.9 - GO Enrichment Analysis with Barplot and Lollipop Plot(GO_Enrichment_Barplot_and_Lollipop)
- This R script performs Gene Ontology (GO) enrichment analysis based on differential gene expression results, and visualizes the top enriched - GO Biological Process terms using a barplot with color gradient and a lollipop-style annotated chart.
- Preprocess gene list: Remove Ensembl version suffix
- 
### 7.10 - Enrichment Bar Plot
- This script generates a bar plot to visualize the top 20 enriched pathways or terms from GO or KEGG enrichment analysis. Each bar represents the number of genes involved in a given pathway, with color intensity indicating statistical significance (−log₁₀ p-value).
- 
### 7.10 -  Volcano Plot + Enrichment Barplot Visualization README
- This script generates a volcano plot and a bidirectional enrichment barplot for visualizing differential gene expression and functional enrichment results. Suitable for publication or report export (e.g., as PDF).

### 7.11 - Synchronized Gene Pathway Plot in R
This R script generates a visualization that simultaneously displays pathways and associated marker or target genes in a coherent graphical format. It is designed for biological interpretation of gene enrichment or pathway analysis, enabling clear, intuitive insights into pathway-specific gene distributions.

### 7.11 - Sankey-Bubble Hybrid Plot in R
This R code creates a combined Sankey and bubble plot for enriched terms or pathway data. It visualizes hierarchical relationships and significance or gene counts together, providing a visually engaging summary of biological enrichment results.

### 7.11 - Excel Table Processing Script in Python (`library.xlsx`)
This Python script automates the cleaning and completion of large `.xlsx` datasets (thousands of rows). It installs required libraries, reads the file, checks column 135 for availability, fills in missing values in column 2 by sequential rules, auto-fills column 46 using data from column 5, removes duplicate rows based on column 1, fills NaN if columns 4–6 are all empty, and finally saves the processed file as `final_library.xlsx`.
