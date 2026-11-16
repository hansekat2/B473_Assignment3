# load libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
print(sys.executable)

# PART 1A: load data
expression_data = pd.read_excel(r"C:\Users\kkbri\OneDrive\B473_Python\Gene_Expression_Data.xlsx")
gene_information = pd.read_csv(r"C:\Users\kkbri\OneDrive\B473_Python\Gene_Information.csv")
sample_information = pd.read_csv(r"C:\Users\kkbri\OneDrive\B473_Python\Sample_Information.tsv", sep="\t")

# show first rows of each dataset
print(expression_data.head())
print(gene_information.head())
print(sample_information.head())
print(sample_information.columns) # shows column names

## PART 1B: Change the sample names from the “Gene_Expression_Data.xlsx”, 
# based upon the phenotype presented in “Sample_Information.tsv” 
# (you can keep unique column names by adding a suffix like "_1", "_2" etc) 
sample_names = []
tumor_count = 0
normal_count = 0

# loop through all sample info columns except the first (Probe_ID)
for col in expression_data.columns[1:]: 
    # finds whether sample is labeled as tumor or normal
    phenotype = sample_information.loc[col, "group"]

    # assigns new name based on phenotype
    if phenotype == "tumor":
        tumor_count +=1
        sample_names.append(f"tumor_{tumor_count}")
    else:
        normal_count += 1
        sample_names.append(f"normal_{normal_count}")
# replaces existing sample names but Probe_ID remains unchanged as the first column
expression_data.columns = ["Probe_ID"] + sample_names
print(expression_data.columns.tolist()) # print list of renamed columns

## PART 1C: Split the merged data from part b, into to 2 parts, 
# based upon their labeled phenotype (ie. tumor or normal)
# creates dataframe of tumor columns
tumor_data = expression_data[
    ["Probe_ID"] +
    [col for col in expression_data.columns 
     if col.startswith("tumor_")]
]

# creates dataframe of normal columns
normal_data = expression_data[
    ["Probe_ID"] +
    [col for col in expression_data.columns 
     if col.startswith("normal_")]
]

# shows rows of each subset
print(tumor_data.head())
print(normal_data.head()) 

## PART 1D: Compute the average expression for each probe from the 2 data sets from part c
# gets lists of tumors and normal sample columns
tumor_cols = [
    col for col in expression_data.columns 
    if col.startswith("tumor_")
]

normal_cols = [
    col for col in expression_data.columns 
    if col.startswith("normal_")
]

# calculate mean expression across all samples
expression_data["tumor_mean"] = expression_data[tumor_cols].mean(axis=1)
expression_data["normal_mean"] = expression_data[normal_cols].mean(axis=1)
print(expression_data[["Probe_ID", "tumor_mean", "normal_mean"]].head())

## PART 1E: Determine the fold change for each probe between the two groups 
# ((Tumour – Control) / Control)
expression_data["fold_change"] = (
    expression_data["tumor_mean"] -
    expression_data["normal_mean"]
) / expression_data["normal_mean"]
print(expression_data[["Probe_ID", "tumor_mean", "normal_mean", "fold_change"]].head())

## PART 1F: Use the data from part e and “Gene_Information.csv” to identify all genes fold
#  change magnitude (absolute value) was greater than 5
expression_data["abs_fold_change"] = expression_data["fold_change"].abs() # calculate absolute value of fold change

# merge expression_data w/ gene_information.csv
merged_data = pd.merge(
    expression_data,
    gene_information,
    on="Probe_ID" # uses Probe_ID as matching key
)
print(merged_data.head())

# select genes w/ abs fold change >5
fold_change_genes = merged_data[
    merged_data["abs_fold_change"] > 5
]
print(fold_change_genes)

## PART 1G: Add a column to the result of part f to include if the gene was higher
#  expressed in “Normal” or “Tumor” samples.
# determine which group has higher avg expression level
fold_change_genes["higher_expression"] = np.where(
    fold_change_genes["fold_change"] > 0,
    "Tumor",
    "Normal"
)
print(fold_change_genes[["Probe_ID", "fold_change", "higher_expression"]])

## PART 2A: Perform exploratory data analysis on the genes from part 1g
# show summary statistics
print(fold_change_genes["fold_change"].describe())

# count genes higher in each group
print(fold_change_genes["higher_expression"].value_counts())

## PART 2B: Create a histogram showing the distribution of the number of differentially
# expressed genes (DEGs) by chromosome.
# count number of DEGs on each chromosome
chromosome_counts = (
    fold_change_genes["Chromosome"]
    .value_counts()
)
print(chromosome_counts)

# create histogram
plt.figure(figsize=(8,5))
plt.hist(
    chromosome_counts,
    bins=10,
    edgecolor="black"
)

plt.title("Distribution of DEGs by chromosome")
plt.xlabel("Number of DEGs per chromosome")
plt.ylabel("Frequency of chromosomes")
plt.tight_layout()
plt.show()

## PART 2C: Make another histogram showing the distribution of DEGs by chromosome segregated
#  by sample type (Normal or Tumor)
# count DEGs per chromosome for genes higher in tumor
tumor_counts = (
    fold_change_genes[fold_change_genes["higher_expression"] == "Tumor"]["Chromosome"]
    .value_counts()
)

# count DEGs per chromosome for genes higher in normal
normal_counts = (
    fold_change_genes[fold_change_genes["higher_expression"] == "Normal"]["Chromosome"]
    .value_counts()
)

# create histogram
fig, axes = plt.subplots(1, 2, figsize=(12,5))
axes[0].hist(
    tumor_counts,
    bins=10,
    edgecolor="black",
    color="red"
)

axes[0].set_title("Tumor DEGs")
axes[0].set_xlabel("DEGs per chromosome")
axes[0].set_ylabel("Frequency")

axes[1].hist(
    normal_counts,
    bins=10,
    edgecolor="black",
    color="blue"
)

axes[1].set_title("Normal DEGs")
axes[1].set_xlabel("DEGs per chromosome")
axes[1].set_ylabel("Frequency")

plt.tight_layout()
plt.show()
## note: All genes meeting the threshold |fold change| > 5 were upregulated (higher) in Tumor samples. 
# No genes met the threshold while being more highly expressed in Normal samples; 
# therefore, the histogram for Normal DEGs contained no observations.

## PART 2D: Create a bar chart showing the percentages of the DEGs that are upregulated
#  (higher) in Tumor samples and down regulated (lower) in Tumor samples.
# count unregulated and down regulated genes
deg_counts = fold_change_genes["higher_expression"].value_counts()

# convert counts to percentages
deg_percentages = (
    deg_counts / deg_counts.sum()
) * 100
print(deg_percentages)

# create bar chart
plt.figure(figsize=(6, 4))
plt.bar(
    deg_percentages.index,
    deg_percentages.values,
    color=["red","blue"]
)

plt.title("Percentage of DEGs by expression direction")
plt.xlabel("Expression category")
plt.ylabel("Percentage")

# add percentage labels above bars
for i, value in enumerate(deg_percentages.values):
    plt.text(
        i,
        value + 1,
        f"{value:.1f}%",
        ha="center"
    )

plt.tight_layout()
plt.show()

## PART 2E: Use the raw data from part 1b to create a heatmap visualizing gene
#  expression by sample.
# get Probe_IDs for significant genes from part 1G
deg_ids = fold_change_genes["Probe_ID"]

# select only genes from the raw expression data
heatmap_data = expression_data[expression_data["Probe_ID"].isin(deg_ids)]

# use Probe_ID as the row labels
heatmap_data = heatmap_data.set_index("Probe_ID")

# select only sample columns
sample_cols = [
    col for col in heatmap_data.columns
    if col.startswith("tumor_") or col.startswith("normal_")
]

# extract expression values only
heatmap_values = heatmap_data[sample_cols]

# standardize each gene(row), mean=0, std=1
heatmap_scaled = (
    heatmap_values.sub(heatmap_values.mean(axis=1), axis=0)
    .div(heatmap_values.std(axis=1), axis=0)
)

# create heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(
    heatmap_scaled,
    cmap="coolwarm",
    center=0
)

plt.title("Scaled heatmap of gene expression by sample")
plt.xlabel("Samples")
plt.ylabel("Genes")
plt.tight_layout()
plt.show()

## PART 2F: Use the same data from the previous part to create a clustermap visualizing
#  gene expression by sample.
# extract DEG expression data from the original dataset
clustermap_data = expression_data[expression_data["Probe_ID"].isin(deg_ids)]

# use Probe_ID as row labels
clustermap_data = clustermap_data.set_index("Probe_ID")

# keep only tumor and normal sample columns
sample_cols = [
    col for col in clustermap_data.columns
    if col.startswith("tumor_") or col.startswith("normal_")
]

# expression values only
expression_values = clustermap_data[sample_cols]

# standardize each gene (row)
expression_scaled = (
    expression_values.sub(expression_values.mean(axis=1), axis=0)
    .div(expression_values.std(axis=1), axis=0)
)

# create clustermap
sns.clustermap(
    expression_scaled,
    cmap="coolwarm",
    center=0,
    figsize=(12,10)
)

plt.suptitle("Clustermap of gene expression by sample", y=1.02)
plt.show()

## PART 2G: 45 DEGs were identified using an absolute fold change threshold greater than 5. 
# All of the identified DEGs were upregulated (higher) in Tumor samples, meaning none were
# upregulated in the Normal samples. You can see this in the bar chart which shows that 100% of the
# significant DEGs were associated with higher expression in Tumor samples. The scaled
# heatmap of gene expression revealed the differences between Tumor and Normal samples.
# Many genes showed higher expression in Tumor samples (warm colors), while the Normal samples
# showed lower expression levels.