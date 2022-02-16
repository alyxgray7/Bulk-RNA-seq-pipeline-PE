### Set up
##########
# command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

# help function
help <- function() {
  cat("deseq2_group.R : Performs log-ratio test (LRT) for differential expression between sample groups.
      ")
  cat("\n")
  cat("Usage : \n")
  cat("--countsFile       : [ required ] Filtered counts table.
      \n")
  cat("--outDir           : [ required ] Directory to save outputs to.
    \n")
  cat("--metaFile         : [ required ] Metadata file as specified by config['omic_meta_data']
  \n")
  cat("--sampleID         : [ required ] Column name in config['omic_meta_data'] that identifies sample names in the counts table.
  \n")
  cat("--linear_model     : [ required ] Column name that contains sample group identifiers to build contrast. Also specified by config['linear_model']
  \n")
  cat("--plotCols         : [ required ] Column names in metadata file to annotate pheatmaps by. As specified by config['meta_columns_to_plot']
  \n")
  cat("--LRT              : [ required ] Subset of samples to run DESeq2 on (> pairwise). As specified by config['diffexp']['LRT'].
  \n")
  cat("--pca_labels       : [ required ] Sample identifiers to label PCA plots by. As specified by config['pca']['labels']
  \n")
  cat("\n")
  q()
}

# save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  countsFile <- sub('--countsFile=', '', args[grep('--countsFile=', args)])
  outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
  metaFile <- sub('--metaFile=', '', args[grep('--metaFile=', args)])
  sampleID <- sub('--sampleID=', '', args[grep('--sampleID=', args)])
  linear_model <- sub('--linear_model=', '', args[grep('--linear_model=', args)])
  plotCols <- sub('--plotCols=', '', args[grep('--plotCols=', args)])
  LRT <- sub("--LRT=", '', args[grep('--LRT=', args)])
  pca_labels <- sub("--pca_labels=", '', args[grep("--pca_labels=", args)])
}

# make arguments R-readable
plotCols <- "Lane,Group,Sex,Age_at_collection,Novogene_RIN,TM_class RINcat"
plotCols <- c(strsplit(plotCols, split = ",", fixed = TRUE)[[1]])

LRT <- "1_Case 2_InSitu 3_Control 4_ScreenNegative"
LRT <- c(strsplit(LRT, split = " ", fixed = TRUE)[[1]])

pca_labels <- "Lane Group Sex Age_at_collection Novogene_RIN TM_class RINcat"
pca_labels <- c(strsplit(pca_labels, split = " ", fixed = TRUE)[[1]])

# save std in/out
io <- list(
  countsFile = countsFile,
  outDir = outDir,
  metaFile = metaFile,
  sampleID = sampleID,
  Type = linear_model,
  group = group,
  plot_cols = plotCols,
  labels = pca_labels
)
io

# cat(sprintf(c('Working directory',getwd(), '\n')))
# cat(sprintf('Setting parameters', '\n'))

# debug on exa
io <- list(
  countsFile = "data/platelet_full-cohort_genecounts.filt.txt"
  , outDir = "results/diffexp/group"
  , metaFile = "data/md_merged_noI9.tsv"
  , sampleID = "rnaSampleID"
  , Type = "Group"
  , group = LRT
  , plot_cols = plotCols
  , labels = pca_labels
  , subset_cols = plotCols
)
io

### Inputs
##########
# labels <- snakemake@params[['pca_labels']]
# labels <- "Lane Group Sex Age_at_collection Novogene_RIN TM_class RINcat"
# cat(sprintf(c('PCA Labels: ',labels, "\n")))

# counts <- snakemake@input[['counts']]
# cat(sprintf(c('Counts table: ', counts, '\n')))

# metadata <- snakemake@params[['samples']]
# cat(sprintf(c('Metadata: ', metadata, '\n')))
# metadata <- "data/md_merged_noI9.tsv"

# sampleID <- snakemake@params[['sample_id']]
# cat(sprintf(c('Sample ID: ', sampleID, '\n')))
# sampleID <- "rnaSampleID"

# Type <- snakemake@params[['linear_model']]
# cat(sprintf(c('Linear Model: ', Type, '\n')))
# Type <- "Group"

# group <- snakemake@params[['LRT']]
# cat(sprintf(c('Subsetted group: ', group, '\n')))
# group <- "1_Case 2_InSitu 3_Control 4_ScreenNegative"

# # plot_cols <- snakemake@config[['meta_columns_to_plot']]
# plot_cols <- "Lane Group Sex Age_at_collection Novogene_RIN TM_class RINcat"

# subset_cols = names(plot_cols)
# subset_cols <- "Lane Group Sex Age_at_collection Novogene_RIN TM_class RINcat"


### Outputs
###########
# Dir <- "results/diffexp/group/"

# pca_plot <- snakemake@output[['pca']]
# pca_plot <- "results/diffexp/group/LRT_pca.pdf"
pca_plot <- paste0(io$outDir, "/", "LRT_pca.pdf")
cat(sprintf(c('PCA plot: ', pca_plot, "\n")))

# sd_mean_plot <- snakemake@output[['sd_mean_plot']]
sd_mean_plot <- paste0(io$outDir, "/", "LRT_sd_mean_plot.pdf")
cat(sprintf(c('SD Mean plot: ',sd_mean_plot,'\n')))
# sd_mean_plot <- "results/diffexp/group/LRT_sd_mean_plot.pdf"

# distance_plot <- snakemake@output[['distance_plot']]
distance_plot <- paste0(io$outDir, "/", "LRT_distance_plot.pdf")
cat(sprintf(c('Distance plot: ',distance_plot,'\n')))
# distance_plot <- "results/diffexp/group/LRT_distance_plot.pdf"

# heatmap_plot <- snakemake@output[['heatmap_plot']]
heatmap_plot <- paste0(io$outDir, "/", "LRT_heatmap_plot.pdf")
cat(sprintf(c('Heatmap Plot: ', heatmap_plot, '\n')))
# heatmap_plot <- "results/diffexp/group/LRT_heatmap_plot.pdf"

# rds_out <- snakemake@output[['rds']]
rds_out <- paste0(io$outDir, "/", "LRT_all.rds")
cat(sprintf(c('RDS Output: ', rds_out, '\n')))
# rds_out <- "results/diffexp/group/LRT_all.rds"

# rld_out <- snakemake@output[['rld_out']]
rld_out <- paste0(io$outDir, "/", "LRT_rlog_dds.rds")
cat(sprintf(c('RLD Output: ', rld_out, '\n')))
# rld_out <- "results/diffexp/group/LRT_rlog_dds.rds"


# color palette
# colors <- snakemake@params[['colors']]
colors <- "NA"
# discrete <- snakemake@params[['discrete']]
discrete <- "NA"

# libraries
library("DESeq2")
library("ggplot2")
library("pheatmap")
library("dplyr")
library("vsn")
library("RColorBrewer")
library("genefilter")

# function to grab the ggplot2 colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)
md <- read.delim(io$metaFile, sep = "\t", stringsAsFactors = FALSE)
# md <- md[order(md[io$sampleID]),]
md <- md[order(md[[io$sampleID]]), ]
head(md)

# Read in counts table
# cts <- read.table(counts, header=TRUE, row.names=1, sep="\t", check.names=F)
cts <- read.table(io$countsFile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cts <- cts[, order(colnames(cts))]
dim(cts)

# Put sample IDs as rownames of metadata
rownames(md) <- md[[io$sampleID]]
md[[io$sampleID]] <- NULL

# Ensure that we subset md to have exactly the same samples as in the counts table
cts <- cts[, colnames(cts) %in% rownames(md)]
md <- md[colnames(cts),]
dim(md)
head(md)

# Check
stopifnot(rownames(md)==colnames(cts))

# Define colours based on number of Conditions
if(colors[[1]] != 'NA' & discrete[[1]] == 'NA') {
    if (brewer.pal.info[colors[[1]],]$maxcolors >= length(unique(md[[io$Type]]))) {
        pal <- brewer.pal(length(unique(md[[io$Type]])), name = colors[[1]])
    } 
} else if(discrete[[1]] != 'NA' & length(discrete) == length(unique(md[[io$Type]]))){
        pal <- unlist(discrete)
} else {
        pal <- gg_color_hue(length(unique(md[[io$Type]])))
}

print(as.formula(paste('~', io$Type)))

# Create dds object from counts data and correct columns
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = md,
                              design = as.formula(paste('~', io$Type)))

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]
dim(dds)

# Likelihood Ratio test to look at differential expression across ALL types, and not just pairs of types (contrast)
dds.lrt <- DESeq(dds, test= "LRT", reduced=~1)
res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering = FALSE)
head(res.lrt)

# Obtain normalized counts
rld <- rlog(dds.lrt, blind = FALSE)

# Pairwise PCA Plot
pdf(pca_plot)
plotPCA(rld, intgroup = io$Type)
dev.off()

# Pairwise PCA Plot with more than one PCA parameter
if (length(io$labels) > 1) {
  pairs <- combn(io$labels, 2)
  for (j in 1:ncol(pairs)) {
    j <- 1
    paired <- paste0(pairs[,j][1], "-vs-", pairs[,j][2])
    # pca_plot2 <- paste0(pairs[,j][1], "-vs-", pairs[,j][2], "_pcaPlot.pdf")
    pca_plot2 <- paste0(io$outDir, "/", paired, "_pcaPlot.pdf")
    pcaData <- plotPCA(rld, intgroup = c(pairs[,j][1], pairs[,j][2]), returnData = TRUE)
    pdf(pca_plot2, 5, 5)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes_string("PC1", "PC2", color = pairs[,j][1], shape = pairs[,j][2])) +
      geom_point(size = 3) +
      ggtitle(paired) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed()
    dev.off()
  }
  # # pca_plot2 <- sub("$", "twoDimensional_pca_plot.pdf", paste0(io$outDir, "/"))
  # # pcaData <- plotPCA(rld, intgroup = c(labels[[1]], labels[[2]]), returnData = TRUE)
  # pdf(pca_plot2, 5, 5)
  # percentVar <- round(100 * attr(pcaData, "percentVar"))
  # ggplot(pcaData, aes_string("PC1", "PC2", color = labels[[1]], shape = labels[[2]])) +
  #   geom_point(size = 3) +
  #   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  #   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #   coord_fixed()
  # dev.off()
}

# SD mean plot
pdf(sd_mean_plot)
meanSdPlot(assay(rld))
dev.off()

# Heatmap of distances
pdf(distance_plot)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- colnames(rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, fontsize=5, scale="row",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# Heatmap across all samples
# List top 50 genes for group comparisons
topGenes <- head(order(res.lrt$padj), 50)

# Extract topGenes from rld object
plot <- assay(rld)[topGenes,] #for 2+ types

# Generate data frame with samples as the rownames and single colData as the first row
# Default when we subset creates an incompatible dataframe so this is a check
df <- as.data.frame(colData(rld))
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot[[sampleID]]
  annot[[sampleID]] <- NULL
} else {
  annot <- df[,subset_cols]
}

pdf(heatmap_plot)
pheatmap(assay(rld)[topGenes,], cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes across all samples"))
dev.off()

saveRDS(dds, file=rds_out)
saveRDS(rld, file=rld_out)

group <- as.vector(group)

# If LRT group has been specified, run the analysis for that group
if (length(group)>0) {
  md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)
  md <- md[order(md[sampleID]),]
  cts <- read.table(counts, header=TRUE, row.names=1, sep="\t")
  cts <- cts[,order(colnames(cts))]
  md <- md[md[[Type]] %in% group,]
  rownames(md) <- md[[sampleID]]
  md[[sampleID]] <- NULL
  keep <- colnames(cts)[colnames(cts) %in% rownames(md)]
  cts <- cts[, keep]
  dim(cts)
  md <- md[colnames(cts),]
  dim(md)

  dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=md,
                              design= as.formula(paste('~',Type)))
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  dds.lrt <- DESeq(dds, test="LRT", reduced=~1)
  res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering=FALSE)
  rld <- rlog(dds.lrt, blind=FALSE)
    
  # Pairwise PCA Plot
  pdf(sub("$", "subsetted_pca_plot.pdf", Dir), 5, 5)
  plotPCA(rld, intgroup=labels[[1]])
  dev.off()
  # Pairwise PCA Plot with more than one PCA parameter
  if (length(labels)>1) {
    pcaData <- plotPCA(rld, intgroup=c(labels[[1]], labels[[2]]), returnData=TRUE)
    pdf(sub("$", "subsetted_twoDimensional_pca_plot.pdf", Dir), 5, 5)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes_string("PC1", "PC2", color=labels[[1]], shape=labels[[2]])) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed()
    dev.off()
  }
  
  # Heatmap
  topGenes <- head(order(res.lrt$padj), 50)
  # Extract topGenes from rld object
  plot <- assay(rld)[topGenes,] #for 2+ types

  df <- as.data.frame(colData(rld))
  if (length(subset_cols)==1) {
    annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
    names(annot) <- c("SampleID", subset_cols[1])
    rownames(annot) <- annot[[sampleID]]
    annot[[sampleID]] <- NULL
  } else {
    annot <- df[,subset_cols]
  }

  pdf(sub("$", "subsetted_heatmap.pdf", Dir), 5, 5)
  pheatmap(assay(rld)[topGenes,], cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes across selected samples"))
  dev.off()
}
