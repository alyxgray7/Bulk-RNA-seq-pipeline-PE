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
  cat("--colors           : [ required ] Preferred colors for ggplot. As specified by config['colors']['rcolorbrewer']
  \n")
  cat("--discrete         : [ required ] Preferred colors for discrete variables. As specified by config['colors']['discrete']
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
  colors <- sub("--colors=", '', args[grep("--colors=", args)])
  discrete <- sub("--discrete=", '', args[grep("--discrete=", args)])
}

# # for debugging on exa
# plotCols <- "lane,group,sex,age,TM_class,RINcat,diabetes"
# LRT <- "1_Case 2_InSitu 3_Control 4_ScreenNegative"
# pca_labels <- "group sex RINcat diabetes"

# make arguments R-readable
plotCols <- c(strsplit(plotCols, split = ",", fixed = TRUE)[[1]])
LRT <- c(strsplit(LRT, split = " ", fixed = TRUE)[[1]])
pca_labels <- c(strsplit(pca_labels, split = " ", fixed = TRUE)[[1]])

# save std in/out
io <- list(
  countsFile = countsFile,
  outDir = outDir,
  metaFile = metaFile,
  sampleID = sampleID,
  Type = linear_model,
  group = LRT,
  plot_cols = plotCols,
  labels = pca_labels,
  subset_cols = plotCols,
  colors = colors,
  discrete = discrete
)
io

# cat(sprintf(c('Working directory',getwd(), '\n')))
# cat(sprintf('Setting parameters', '\n'))

# # debug on exa
# io <- list(
#   countsFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_02162022/data/platelet_full-cohort_rmdp_genecounts.filt.txt"
#   , outDir = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_02162022/results/diffexp/group"
#   , metaFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_02162022/data/pltRNAseq_metadata_02162022.tsv"
#   , sampleID = "rnaSampleID"
#   , Type = "group"
#   , group = LRT
#   , plot_cols = plotCols
#   , labels = pca_labels
#   , subset_cols = plotCols
#   , colors = "NA"
#   , discrete = "NA"
# )
# io

# create outdir if needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}


### Outputs
###########
# pca_plot <- snakemake@output[['pca']]
pca_plot <- paste0(io$outDir, "/", "LRT_pca.pdf")
cat(sprintf(c('PCA plot: ', pca_plot, "\n")))

# sd_mean_plot <- snakemake@output[['sd_mean_plot']]
sd_mean_plot <- paste0(io$outDir, "/", "LRT_sd_mean_plot.pdf")
cat(sprintf(c('SD Mean plot: ',sd_mean_plot,'\n')))

# distance_plot <- snakemake@output[['distance_plot']]
distance_plot <- paste0(io$outDir, "/", "LRT_distance_plot.pdf")
cat(sprintf(c('Distance plot: ',distance_plot,'\n')))

# heatmap_plot <- snakemake@output[['heatmap_plot']]
heatmap_plot <- paste0(io$outDir, "/", "LRT_heatmap_plot.pdf")
cat(sprintf(c('Heatmap Plot: ', heatmap_plot, '\n')))

# rds_out <- snakemake@output[['rds']]
rds_out <- paste0(io$outDir, "/", "LRT_all.rds")
cat(sprintf(c('RDS Output: ', rds_out, '\n')))

# rld_out <- snakemake@output[['rld_out']]
rld_out <- paste0(io$outDir, "/", "LRT_rlog_dds.rds")
cat(sprintf(c('RLD Output: ', rld_out, '\n')))

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


### Format data
###############
# metadata
md <- read.delim(io$metaFile, sep = "\t", stringsAsFactors = FALSE)
md <- md[order(md[[io$sampleID]]), ]
head(md)

# Put sample IDs as rownames of metadata
rownames(md) <- md[[io$sampleID]]
md[[io$sampleID]] <- NULL

# counts table
cts <- read.table(io$countsFile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cts <- cts[, order(colnames(cts))]
dim(cts)

# subset counts to have exactly the same samples as in the md table
cts <- cts[, colnames(cts) %in% rownames(md)]
md <- md[colnames(cts),]
dim(md)
head(md)

# check
stopifnot(rownames(md) == colnames(cts))

# define colours based on number of Conditions
if(io$colors[[1]] != 'NA' & io$discrete[[1]] == 'NA') {
    if (brewer.pal.info[io$colors[[1]],]$maxcolors >= length(unique(md[[io$Type]]))) {
        pal <- brewer.pal(length(unique(md[[io$Type]])), name = io$colors[[1]])
    } 
} else if(io$discrete[[1]] != 'NA' & length(io$discrete) == length(unique(md[[io$Type]]))){
        pal <- unlist(io$discrete)
} else {
        pal <- gg_color_hue(length(unique(md[[io$Type]])))
}


### Run LRT on all
##################
# create dds object from counts data and correct columns
print(as.formula(paste('~', io$Type)))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = md,
                              design = as.formula(paste('~', io$Type)))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]
dim(dds)

# run likelihood ratio test
# looks at differential expression across ALL types, and not just pairs of types (contrast)
dds.lrt <- DESeq(dds, test = "LRT", reduced=~1)
res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering = FALSE)
head(res.lrt)

# save LRT results
res.lrt.df <- as.data.frame(res.lrt)
res.lrt.df[["Gene"]] <- rownames(res.lrt.df)
res.lrt.df <- res.lrt.df[, c(ncol(res.lrt.df), 1:(ncol(res.lrt.df)-1))]
write.table(res.lrt.df, paste0(io$outDir, "/LRT_results.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# obtain rlog normalized counts
rld <- rlog(dds.lrt, blind = FALSE)
write.table(assay(rld), paste0(io$outDir, "/rlog_counts.tsv"), col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


### Summary plots
#################
# pairwise PCA plot
pdf(pca_plot)
plotPCA(rld, intgroup = io$Type)
dev.off()

# pairwise PCA plot with sample names
pcaData <- plotPCA(rld, intgroup = io$Type, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = as.factor(group), label = name)) +
  geom_text() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
ggsave(paste0(io$outDir, "/PCA_text.pdf"), device = "pdf")
  
# pairwise PCA plots with more than one PCA parameter
if (length(io$labels) > 2) {
  pairs <- combn(io$labels, 2)
  for (j in 1:ncol(pairs)) {
    paired <- paste0(pairs[,j][1], "-vs-", pairs[,j][2])
    pcaData <- plotPCA(rld, intgroup = c(pairs[,j][1], pairs[,j][2]), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes_string("PC1", "PC2", color = pairs[,j][2], shape = pairs[,j][1])) +
      geom_point(size = 3) +
      ggtitle(paired) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed()
    ggsave(paste0(io$outDir, "/", paired, "_pcaPlot.pdf"), device = "pdf")
  }
}

if (length(io$labels) == 2){
  pca_plot2 <- sub("$", "twoDimensional_pca_plot.pdf", paste0(io$outDir, "/"))
  pcaData <- plotPCA(rld, intgroup = c(io$labels[[1]], io$labels[[2]]), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes_string("PC1", "PC2", color = io$labels[[1]], shape = io$labels[[2]])) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  ggsave(pca_plot2, device = "pdf")
}

# SD mean plot
pdf(sd_mean_plot)
meanSdPlot(assay(rld))
dev.off()

# sample distances heatmap
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

# heatmap across all samples
# list top 50 genes for group comparisons
topGenes <- head(order(res.lrt$padj), 50)

# extract topGenes from rld object
plot <- assay(rld)[topGenes,] # for 2+ types

# generate data frame with samples as the rownames and single colData as the first row
# default when we subset creates an incompatible dataframe so this is a check
df <- as.data.frame(colData(rld))
if (length(io$subset_cols) == 1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[io$subset_cols[1]]])))
  names(annot) <- c("SampleID", io$subset_cols[1])
  rownames(annot) <- annot[[io$sampleID]]
  annot[[io$sampleID]] <- NULL
} else {
  annot <- df[,io$subset_cols]
}

# annotated clustered heatmap of top DE genes
pdf(heatmap_plot)
pheatmap(assay(rld)[topGenes,], cluster_rows=T, scale="row", fontsize=6,fontsize_row=6, fontsize_col=6, show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes across all samples"))
dev.off()

# save Rdata objects
saveRDS(dds, file=rds_out)
saveRDS(rld, file=rld_out)


### Run LRT subset
##################
# run only for groups included
group <- as.vector(io$group)
if (length(group) > 0) {

  # reload metadata and counts
  md <- read.delim(io$metaFile, sep = "\t", stringsAsFactors = FALSE)
  md <- md[order(md[[io$sampleID]]), ]
  cts <- read.table(io$countsFile, header = TRUE, row.names = 1, sep = "\t")
  cts <- cts[ , order(colnames(cts))]
  md <- md[md[[io$Type]] %in% group, ]
  rownames(md) <- md[[io$sampleID]]
  md[[io$sampleID]] <- NULL
  keep <- colnames(cts)[colnames(cts) %in% rownames(md)]
  cts <- cts[, keep]
  dim(cts)
  md <- md[colnames(cts),]
  dim(md)

  # create deseq2 object
  dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = md,
                              design = as.formula(paste('~', io$Type)))
  
  # remove uninformative rows
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  
  # run LRT
  dds.lrt <- DESeq(dds, test = "LRT", reduced=~1)
  res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering=FALSE)
  
  # save LRT on subsetted samples
  res.lrt.df <- as.data.frame(res.lrt)
  res.lrt.df[["Gene"]] <- rownames(res.lrt.df)
  res.lrt.df <- res.lrt.df[, c(ncol(res.lrt.df), 1:(ncol(res.lrt.df)-1))]
  write.table(res.lrt.df, paste0(io$outDir, "/subsetted_LRT_results.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

  # obtain rlog normalized counts
  rld <- rlog(dds.lrt, blind = FALSE)
  write.table(assay(rld), paste0(io$outDir, "/subsetted_rlog_counts.tsv"), col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
    

  ### Summary plots
  #################
  # pairwise PCA plot
  pdf(paste0(io$outDir, "/subsetted_PCA.pdf"))
  plotPCA(rld, intgroup = io$Type)
  dev.off()

  # pairwise PCA plot with sample labels
  pcaData <- plotPCA(rld, intgroup = io$Type, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(x = PC1, y = PC2, color = as.factor(group), label = name)) +
    geom_text() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  ggsave(paste0(io$outDir, "/subsetted_PCA_text.pdf"), device = "pdf")

  # pairwise PCA plot with sample names
  pcaData <- plotPCA(rld, intgroup = io$Type, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(x = PC1, y = PC2, color = as.factor(group), label = name)) +
    geom_text() +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  ggsave(paste0(io$outDir, "/subsetted_PCA_text.pdf"), device = "pdf")
    
  # pairwise PCA plots with more than one PCA parameter
  if (length(io$labels) > 2) {
    pairs <- combn(io$labels, 2)
    for (j in 1:ncol(pairs)) {
      paired <- paste0(pairs[,j][1], "-vs-", pairs[,j][2])
      pcaData <- plotPCA(rld, intgroup = c(pairs[,j][1], pairs[,j][2]), returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      ggplot(pcaData, aes_string("PC1", "PC2", color = pairs[,j][2], shape = pairs[,j][1])) +
        geom_point(size = 3) +
        ggtitle(paired) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed()
      ggsave(paste0(io$outDir, "/subsetted_", paired, "_pcaPlot.pdf"), device = "pdf")
    }
  }

  if (length(io$labels) == 2){
    pca_plot2 <- sub("$", "subsetted_twoDimensional_pca_plot.pdf", paste0(io$outDir, "/"))
    pcaData <- plotPCA(rld, intgroup = c(io$labels[[1]], io$labels[[2]]), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes_string("PC1", "PC2", color = io$labels[[1]], shape = io$labels[[2]])) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed()
    ggsave(pca_plot2, device = "pdf")
  }
  
  # heatmap
  topGenes <- head(order(res.lrt$padj), 50)

  # extract topGenes from rld object
  plot <- assay(rld)[topGenes,] #for 2+ types

  df <- as.data.frame(colData(rld))
  if (length(io$subset_cols) == 1) {
    annot <- as.data.frame(cbind(rownames(df), paste(df[[io$subset_cols[1]]])))
    names(annot) <- c("SampleID", io$subset_cols[1])
    rownames(annot) <- annot[[io$sampleID]]
    annot[[io$sampleID]] <- NULL
  } else {
    annot <- df[,io$subset_cols]
  }

  # clustered heatmap
  pdf(sub("$", "/subsetted_heatmap.pdf", paste0(io$outDir)), 5, 5)
  pheatmap(assay(rld)[topGenes,], cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6, show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes across selected samples"))
  dev.off()
}
