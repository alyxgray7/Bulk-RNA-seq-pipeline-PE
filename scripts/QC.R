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
  cat("--rld              : [ required ] Rdata object produced from deseq2 and contains the raw counts table.
  \n")
  cat("--rds              : [ required ] Rdata object containing rlog normalized counts from deseq2.")
  cat("--outDir           : [ required ] Directory to save outputs to.
    \n")
  cat("--sampleID         : [ required ] Column name in config['omic_meta_data'] that identifies sample names in the counts table.
  \n")
  cat("--Type             : [ required ] Column name that contains sample group identifiers to build contrast. Also specified by config['linear_model']
  \n")
  cat("--plot_cols         : [ required ] Column names in metadata file to annotate pheatmaps by. As specified by config['meta_columns_to_plot']
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
  rld <- sub("--rld=", '', args[grep("--rld=", args)])
  rds <- sub('--rds=', '', args[grep("--rds=", args)])
  outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
  sampleID <- sub('--sampleID=', '', args[grep('--sampleID=', args)])
  Type <- sub('--Type=', '', args[grep('--Type=', args)])
  plot_cols <- sub('--plot_cols=', '', args[grep('--plot_cols=', args)])
  colors <- sub("--colors=", '', args[grep("--colors=", args)])
  discrete <- sub("--discrete=", '', args[grep("--discrete=", args)])
}

# plot_cols <- "group,diabetes,sex,age,TM_class,RINcat"
plot_cols <- c(strsplit(plot_cols, split = ",", fixed = TRUE)[[1]])

io <- list(
  rld = rld
  , rds = rds
  , outDir = outDir
  , sampleID = sampleID
  , Type = Type
  , plot_cols = plot_cols
  , subset_cols = plot_cols
  , colors = colors
  , discrete = discrete
)

# debugging on exa
io <- list(
  rld = "results/diffexp/group/LRT_rlog_dds.rds"
  , rds = "results/diffexp/group/LRT_all.rds"
  , outDir = "results/diffexp/group_test"
  , sampleID = "rnaSampleID"
  , Type = "group"
  , plot_cols = c(plot_cols)
  , subset_cols = c(plot_cols)
  , colors = "NA"
  , discrete = "NA"
)
io

# libraries
library("DESeq2")
library("reshape2")
library("cowplot")
library("limma")
library("vsn")
library("genefilter")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("pheatmap")
library("hexbin")

# output files
# MDS_out <- snakemake@output[['mds_plot']]
MDS_out <- paste0(io$outDir, "/MDS_plot.pdf")
print(MDS_out)

# MDS_table <- snakemake@output[['mds_table']]
MDS_table <- paste0(io$outDir, "/MDS_table.txt")
print(MDS_table)

# heatmap_out <- snakemake@output[['heatmap_plot']]
heatmap_out <- paste0(io$outDir, "/Heatmap_all_genes.pdf")
print(heatmap_out)

# sd_out <- snakemake@output[['sd_plot']]
sd_out <- paste0(io$outDir, "/stdev_plot.pdf")
print(sd_out)

# normCounts_out <- snakemake@output[['rlogCounts_plot']]
normCounts_out <- paste0(io$outDir, "/rlog_counts_violinPlot.pdf")
print(normCounts_out)

# normCounts_fac <- snakemake@output[['rlogCounts_fac_plot']]
normCounts_fac <- paste0(io$outDir, "/rlog_counts_faceted_violinPlot.pdf")
print(normCounts_fac)

# rawCounts_out <- snakemake@output[['counts_plot']]
rawCounts_out <- paste0(io$outDir, "/counts_violinPlot.pdf")
print(rawCounts_out)

# rawCounts_fac <- snakemake@output[['counts_fac_plot']]
rawCounts_fac <- paste0(io$outDir, "/counts_faceted_violinPlot.pdf")
print(rawCounts_fac)


### Functions
#############
# function to grab the ggplot2 colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# function to save pheatmap plots
save_pheatmap_pdf <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


### Begin QC
############
# load data
rld <- readRDS(io$rld)
dds <- readRDS(io$rds)

# extract raw counts
rawCounts <- counts(dds, normalized=FALSE)
dim(rawCounts)
head(rawCounts)

# extract metadata
md <- as.data.frame(colData(rld))
md$SampleID <- rownames(md)
head(md)

# define colors
if(io$colors[[1]] != 'NA' & io$discrete[[1]] == 'NA'){
    if (brewer.pal.info[io$colors[[1]],]$maxcolors >= length(unique(md[[io$Type]]))) {
        pal <- brewer.pal(length(unique(md[[io$Type]])), name=io$colors[[1]])
    } 
} else if(io$discrete[[1]] != 'NA' & length(io$discrete)==length(unique(md[[io$Type]]))){
        pal <- unlist(io$discrete)
} else {
        pal <- gg_color_hue(length(unique(md[[io$Type]])))
}

# violin plot of raw counts
df1 <- melt(rawCounts) %>%
  dplyr::rename(Gene=Var1) %>%
  dplyr::rename(SampleID=Var2) %>%
  dplyr::rename(counts=value)
iv <- match(df1$SampleID, md$SampleID)
df1$Condition <- paste(md[iv,][[io$Type]])
df1$SampleID <- factor(df1$SampleID, levels=unique(md$SampleID))
head(df1)
dim(df1)

# aesthetic for plots
dodge <- position_dodge(width = 0.6)
theme_update(plot.title = element_text(hjust = 0.5))

# plot!
p1 <- ggplot(data=df1, mapping=aes(x=SampleID, y=counts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_y_log10() +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=6))

# width of pdf to ensure all sampleIDs are visible when exported to pdf
# This was generated with a use case of 16 samples and a width of 7 fitting well, the +8 is to account for the margins
width <- 7/24*(nrow(md)+8)

# raw counts boxplot
pdf(rawCounts_out, width, 5)
print({
  p1
})
dev.off()

# faceted by condition
p2 <- ggplot(data=df1, mapping=aes(x=SampleID, y=counts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_y_log10() +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=4)) +
  facet_wrap(~Condition)

pdf(rawCounts_fac, 2*width, 5)
print({
  plot_grid(p1, p2)
})
dev.off()

# run same analysis for log2-transformed normalized counts
df2 <- melt(assay(rld)) %>%
  dplyr::rename(Gene=Var1) %>%
  dplyr::rename(SampleID=Var2) %>%
  dplyr::rename(normCounts=value)

# add condition information to this dataframe
iv <- match(df2$SampleID, md$SampleID)
df2$Condition <- paste(md[iv,][[io$Type]])
df2$SampleID <- factor(df2$SampleID, levels=unique(md$SampleID))
head(df2)
dim(df2)

p1 <- ggplot(data=df2, mapping=aes(x=SampleID, y=normCounts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=6)) +
  ylab("regularized log expression")

# raw counts boxplot
pdf(normCounts_out, width, 5)
print({
  p1
})
dev.off()

# faceted by condition
p2 <- ggplot(data=df2, mapping=aes(x=SampleID, y=normCounts, fill=Condition)) +
  geom_violin(width=0.7) +
  geom_boxplot(width=0.2, outlier.colour=NA, position = dodge, color="gray28") +
  scale_fill_manual(values=pal) +
  theme(axis.text.x = element_text(hjust=1, angle=45, size=4)) +
  facet_wrap(~Condition) +
  ylab("regularized log expression")

pdf(normCounts_fac, 2*width, 5)
print({
  plot_grid(p1, p2)
})
dev.off()

# Standard deviation vs. mean
ntd <- normTransform(dds)

pdf(sd_out)
meanSdPlot(assay(ntd))
dev.off()

# Generate annotation column for heatmap
if (length(io$subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(md), paste(md[[io$subset_cols[1]]])))
  names(annot) <- c("SampleID", io$subset_cols[1])
  rownames(annot) <- annot$SampleID
  annot$SampleID <- NULL
} else {
  annot <- md[,io$subset_cols]
}
head(annot)
dim(annot)

# remove 0's from counts table
if ( table(assay(rld) == 0)[["TRUE"]] > 0) {
  mat <- assay(rld)
  keep <- apply(mat, 1, function(x) any(x != 0))
  mat <- mat[keep,]
} else {
  mat <- assay(rld)
}

# clustered heatmap
hm <- pheatmap(
  mat, 
  show_rownames=F, 
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "correlation", 
  clustering_method = "average", 
  annotation_col = annot, 
  scale = "row", 
  main="Unsupervised heatmap of all gene counts across samples",
  fontsize_row=4, fontsize_col=6, fontsize=8,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
save_pheatmap_pdf(hm, heatmap_out)

# use plotMA function from limma, then extract data from this variable to plot with ggplot2
p <- plotMDS(assay(rld), top = 1000)
# df <- data.frame(x=p$x, y=p$y, name=names(p$x))
# df <- data.frame(x = p$x, y = p$y, name=attr(p[[5]], "dimnames")[[1]])
df <- data.frame(x = p$x, y = p$y, name = dimnames(p[[5]])[[1]])
iv <- match(df$name, md$SampleID)
df$Condition <- paste(md[iv,][[io$Type]])

pdf(MDS_out)
ggplot(data=df, mapping=aes(x=x,y=y)) +
  geom_point(size=3, colour = "black", show.legend = TRUE) +
  geom_point(aes(color=Condition), size=2.2) +
  scale_colour_manual(values=pal) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  ggtitle("MDS Plot")
dev.off()

# Export the table for MDS
write.table(df, file=MDS_table, sep="\t", quote=F, row.names=FALSE)

