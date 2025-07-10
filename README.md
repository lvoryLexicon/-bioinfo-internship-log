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
### 7.10 - Enrichment Bar Plot
- This script generates a bar plot to visualize the top 20 enriched pathways or terms from GO or KEGG enrichment analysis. Each bar represents the number of genes involved in a given pathway, with color intensity indicating statistical significance (−log₁₀ p-value).
### 7.10 -  Volcano Plot + Enrichment Barplot Visualization README
- This script generates a volcano plot and a bidirectional enrichment barplot for visualizing differential gene expression and functional enrichment results. Suitable for publication or report export (e.g., as PDF).
