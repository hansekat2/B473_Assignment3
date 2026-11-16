# B473_Assignment3
Programmer: Katie Hansen
Course: INFO-B473
Programming Language: Python 3.14
Date Submitted: 11/15/26

# Purpose
The purpose of this script is to analyze gene expression data from tumor and normal tissue samples. The script renames samples based on phenotype, calculates average expression values, determines fold changes between tumor and normal groups, identifies differentially expressed genes (DEGs), performs exploratory data analysis, and generates multiple visualizations including histograms, bar charts, heatmaps, and a clustermap.

# Required Input Files
1. Gene_Expression_Data.xlsx
2. Gene_Information.csv
3. Sample_Information.tsv

# Required Libraries & Software
Python 3.14
- numpy
- pandas
- matplotlib
- seaborn
- scipy

# Instructions for running the script
Download all required input files and put them in the same directory. Install all required packages.

# Output Files & Results
## Console Output:
- Summary statistics for fold changes
- Lists of significant DEGs
- DEG counts by chromosome
- DEG expression category counts
## Figures Displayed
- Distribution of DEGs by chromosome histogram
- Tumor versus Normal DEG chromosome histograms
- Percentage of DEGs by expression direction bar chart
- Scaled heatmap of gene expression by sample
- Clustermap of gene expression by sample
