# Script for plotting the following pre-DESeq2 QC/QA metrics:
#   1. Sum counts / sample
#   2. Biotype distributions
#   3. Gene attributes

# Set up libraries
library(ggplot2)
library(ggpubr)
library(data.table)
library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(scales)

# Set up variables
print("Setting parameters")

counts = snakemake@input[['counts']]
rDist = snakemake@input[['rDist']]

out_plot = snakemake@output[['combo_plot']] # A. sumCount; B. biotype; C. gene attribute
biotype_table = snakemake@output[['biotype_table']]

md <- snakemake@params[['samples']]
anno <- snakemake@param[['anno']]

print("Passes test")