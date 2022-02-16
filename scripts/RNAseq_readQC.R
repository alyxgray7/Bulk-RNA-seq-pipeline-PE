
### Set up
##########
# command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
help <- function(){
  cat("RNAseq_readQC.R :
  
  - Outputs :
  
  1) ggplot2 figure displaying the number of input reads and percent mitochondrial reads by sample. Also includes violin/boxplot distributions by contrast. 
          - file extension {outdir}/{outname}_readSummary.pdf
  2) ggplot2 figure displaying the distribution of reads mapping to the three gene attributes (exon, intron, intergenic). Plotted by sample.
          - file extension {outdir}/{outname}_mappings.pdf
  3) ggplot2 figure displaying the distribution of reads mapping to the top 6 biotypes by sample.
          - file extension {outdir}/{outname}_biotypes.pdf
  4) Pairwise (sample-sample) raw count correlation plot
          - file extension {outdir}/{outname}_pairCorrplot.png
  5) Sample correlation heatmap with heirarchical clustering, annotated by the metadata
          - file extension {outdir}/{outname}_corrHeatmap.png
      ")
  cat("\n")
  cat("Usage : \n")
  cat("--annoFile     : Path to config['filter_anno']                           [ required ]
      \n")
  cat("--metaFile     : Path to config['omic_meta_data']                        [ required ]
      \n")
  cat("--countsFile   : Path to unfiltered counts matrix                        [ required ]
      \n")
  cat("--readDistFile : Path to read_coverage.txt                               [ required ]
                        (output from rule compile_rd)
      \n")
  cat("--contrast     : The column name in your config['omic_meta_data'] file,  [ required ]
                        this is the characteristic you would like to do DE on.
                        Example: diagnosis, geotype, etc. (used to color the plots by.)
      \n")
  cat("--plotCols     : Column names in metadata file to annotate clustering;   [ required ]
                        specified by config['meta_columns_to_plot']
      \n")
  cat("--outdir       : Path to the directory to save the figures in.           [ required ]
      \n")
  cat("--corType      : Type of correlation test to run.                        [ required ]
                        Options are Pearson, Spearman, or both
      \n")
  cat("\n")
  q()
}

# save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  annoFile     <- sub('--annoFile=', '', args[grep('--annoFile=', args)])
  metaFile     <- sub('--metaFile=', '', args[grep('--metaFile=', args)])
  countsFile   <- sub('--countsFile=', '', args[grep('--countsFile=', args)])
  readDistFile <- sub('--readDistFile=', '', args[grep('--readDistFile=', args)])
  contrast     <- sub('--contrast=', '', args[grep('--contrast=', args)])
  plotCols     <- sub('--plotCols=', '', args[grep('--plotCols=', args)])
  outDir       <- sub('--outdir=', '', args[grep('--outdir=', args)])
  corType      <- sub('--corType=', '', args[grep('--corType=', args)])
  sampleID     <- sub('--sampleID=', '', args[grep('--sampleID=', args)])
}

# make plotCols argument R-readable
plotCols <- c(strsplit(plotCols, split = ",", fixed = TRUE)[[1]])

# check input files in logs/.out file
io <- list(
  annoFile = annoFile,
  metaFile = metaFile,
  countsFile = countsFile,
  readDistFile = readDistFile,
  contrast = contrast,
  plotCols = plotCols,
  outDir = outDir,
  corType = corType,
  sampleID = sampleID
)
io

# ### Debugging on exa
# ####################
# io <- list(
#   annoFile = "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_94.biomaRt.geneAnno.Rdata",
#   metaFile = "data/md_merged_noI9.tsv",
#   # countsFile = "data/platelet_full-cohort_counts.txt",
#   # countsFile = "data/platelet_full-cohort_genecounts.txt",
#   countsFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/data/platelet_full-cohort_genecounts_noERCC.txt",
#   # readDistFile = "results/tables/read_coverage.txt",
#   readDistFile = "results/tables/read_coverage_edit.txt",
#   contrast = "Group",
#   plotCols = c("Lane", "Group", "Sex", "Age_at_collection", "Novogene_RIN", "TM_class"),
#   outDir = "results/readQC_htseqCounts_noI9",
#   corType = "Spearman",
#   sampleID = "rnaSampleID"
# )
# io

# create outdir as needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}

# libraries
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(corrplot)
library(pheatmap) 
library(Hmisc)

### Load and organize data
##########################
# annotation file
load(io$annoFile, verbose = TRUE)

# meta data
md  <- read.table(io$metaFile, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
rownames(md) <- md[[io$sampleID]]
head(md)

# organize sampleIDs so the order is the same
md <- md[order(md[, io$contrast], rownames(md)), ]
#rownames(md) <- md$Alias
ord <- rownames(md)
ord

# unfiltered counts table
raw <- read.table(io$countsFile, sep = "\t", header = TRUE, row.names = 1)
raw <- raw[, colnames(raw) %in% rownames(md)]
raw <- raw[, ord]
dim(raw)

# remove uninformative columns
raw <- raw[which(rowSums(raw) >= 1), ]
dim(raw)

# sample read coverages
rDist <- read.table(io$readDistFile, header = TRUE, sep = "\t", row.names = 1)
rDist <- rDist[rownames(rDist) %in% rownames(md), ]
rDist <- rDist[ord, ]
head(rDist)

### Summarize gene counts
#########################
# set up dataframe 
sumCounts.df          <- as.data.frame(apply(raw, 2, sum))
names(sumCounts.df)   <- "sumCounts"
sumCounts.df$SampleID <- rownames(sumCounts.df)
iv                    <- match(rownames(sumCounts.df), rownames(md))
sumCounts.df$contrast <- md[iv, io$contrast]
stopifnot(rownames(sumCounts.df) == rownames(md)[iv])

# save the order
sumCounts.df$SampleID <- factor(sumCounts.df$SampleID, levels = sumCounts.df$SampleID)
# sumCounts.df$SampleID <- factor(sumCounts.df$SampleID, 
#                                 levels = paste(sumCounts.df[order(sumCounts.df$contrast, -sumCounts.df$sumCounts),"SampleID"])) # keeps the correct numerical order for ggplot
sumCounts.df$contrast <- md[iv, io$contrast]

if (io$contrast == "RINcat") {
  sumCounts.df$contrast <- factor(sumCounts.df$contrast, levels = c("low", "med", "high"))
}

# barplot of sumCounts by sample
totalReads.plot <- ggplot(
    sumCounts.df, 
    aes(x = SampleID, y = sumCounts, fill = contrast)) +
    geom_col(color="black") +
    ylab("Total gene counts") +
    xlab("Sample") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scientific) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = "black"))
totalReads.plot
ggsave(filename = paste(io$outDir, "totalReads_barplot.png", sep = "/"), device = "png")

# ordered barplot by group and sumCount
toplot <- sumCounts.df
toplot <- toplot[order(toplot$contrast, toplot$sumCounts), ]
toplot$SampleID <- factor(toplot$SampleID, levels = toplot$SampleID)

# barplot of sumCounts by sample
totalReads.plot <- ggplot(
    toplot, 
    aes(x = SampleID, y = sumCounts, fill = contrast)) +
    geom_col(color="black") +
    ylab("Total gene counts") +
    xlab("Sample") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scientific) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = "black"))
totalReads.plot
ggsave(filename = paste(io$outDir, "totalReads_barplot_ordered.png", sep = "/"), device = "png")

# ordered barplot by RNA batch order and sumCount
# batch <- tstrsplit(toplot$SampleID, split = "_", fixed = TRUE)[[1]]
# batch <- gsub("*[0-9]", "", batch)
# toplot$batch <- batch
# toplot <- toplot[order(toplot$batch, toplot$sumCount), ]
# toplot$SampleID <- factor(toplot$SampleID, levels = toplot$SampleID)
# totalReads.plot <- ggplot(
#     toplot, 
#     aes(x = SampleID, y = sumCounts, fill = contrast)) +
#     geom_col(color="black") +
#     ylab("Total gene counts") +
#     xlab("Sample") +
#     scale_fill_brewer(palette = "Paired") +
#     scale_y_continuous(labels = scientific) +
#     theme_bw()+
#     theme(axis.text.x = element_text(angle = 90, size = 8),
#         axis.text.y = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 10),
#         axis.line = element_line(colour = "black"))
# totalReads.plot
# ggsave(filename = paste(io$outDir, "totalReads_barplot_orderedBatch.png", sep = "/"), device = "png")

# violin plot of gene counts by sample group
totalReads.boxplot <- ggplot(
    sumCounts.df,
    aes(x = contrast, y = sumCounts, fill = contrast)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(color = "black", width = 0.1, fill = "white") +
    xlab("Contrast") +
    ylab("Total gene counts") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scientific) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(colour = "black"))
totalReads.boxplot
ggsave(filename = paste(io$outDir, "totalReads_boxplot.png", sep = "/"), device = "png")


### Summarize mitochondrial counts
##################################
# check how gene names are annotated in counts
annoType <- c()
if (length(grep("ENSG", rownames(raw))) == 0) {
  print("counts are external gene ids")
  annoType <- "external_gene_name"
  mtGenes <- anno[grep("^MT-", anno$external_gene_name), annoType]
  mtSub.mat <- raw[which(rownames(raw) %in% mtGenes), ]
  #mtSub.mat             <- raw[grep("^MT-", anno$external_gene_name), annoType, ]
  #mtSub.mat             <- raw[grep("^MT-", rownames(raw), value = TRUE, ignore.case = TRUE), ]
  #sumCounts.df$mtCounts <- colSums(mtSub.mat)
  #sumCounts.df$mtCounts <- apply(mtSub.mat, 2, sum)
  #sumCounts.df$fracMT   <- sumCounts.df$mtCounts / sumCounts.df$sumCounts
  sumCounts.df$fracMT <- colSums(mtSub.mat) / colSums(raw)
} else {
  print("counts are ensembl ids")
  annoType <- "ensembl_gene_id"
  ids                   <- anno[grep("^MT-", anno$external_gene_name), "ensembl_gene_id"]
  raw$ids               <- sub("\\..*$", '', rownames(raw)) # be sure ens Ids are unique
  mtSub.mat             <- raw[raw$ids %in% ids, ]
  raw$ids               <- NULL
  mtSub.mat$ids         <- NULL
  sumCounts.df$mtCounts <- colSums(mtSub.mat)
  sumCounts.df$fracMT   <- sumCounts.df$mtCounts / sumCounts.df$sumCounts
}
annoType

# barplot of mito read fraction by sample
mtReads.plot <- ggplot(
    sumCounts.df, 
    aes(x = SampleID, y = fracMT, fill = contrast)) +
    geom_col(color="black") +
    ylab("Percent mitochondrial counts") +
    xlab("Sample") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(colour = "black"))
mtReads.plot
ggsave(filename = paste(io$outDir, "mtReads_barplot.png", sep = "/"), device = "png")

# order samples by group and fracMT
toplot <- sumCounts.df
toplot <- toplot[order(toplot$contrast, toplot$fracMT), ]
toplot$SampleID <- factor(toplot$SampleID, levels = toplot$SampleID)
# sumCounts.df <- sumCounts.df[order(sumCounts.df$contrast, sumCounts.df$fracMT), ]
# sumCounts.df$SampleID <- factor(sumCounts.df$SampleID, levels = sumCounts.df$SampleID)

# barplot of mito read fraction by sample
mtReads.plot <- ggplot(
    toplot, 
    aes(x = SampleID, y = fracMT, fill = contrast)) +
    geom_col(color="black") +
    ylab("Percent mitochondrial counts") +
    xlab("Sample") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(colour = "black"))
mtReads.plot
ggsave(filename = paste(io$outDir, "mtReads_barplot_ordered.png", sep = "/"), device = "png")

# violin plot of mito read fraction by sample group
mtReads.boxplot <- ggplot(
    sumCounts.df,
    aes(x = contrast, y = fracMT, fill = contrast)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(color = "black", width = 0.1, fill = "white") +
    xlab("Contrast") +
    ylab("Percent mitochondrial counts") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(colour = "black"))
mtReads.boxplot
ggsave(filename = paste(io$outDir, "mtReads_boxplot.png", sep = "/"), device = "png")

# total read plots arranged
ggarrange(totalReads.plot, totalReads.boxplot, ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "bottom")
ggsave(filename = paste(io$outDir, "totalReads_arranged.png", sep = "/"), device = "png")

# mt read plots arranged
ggarrange(mtReads.plot, mtReads.boxplot, ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = "bottom")
ggsave(filename = paste(io$outDir, "mtReads_arranged.png", sep = "/"), device = "png")


### Summarize RNA biotypes
##########################
# get gene names from counts matrix
if (!annoType == "ensembl_gene_id") {
  fD <- as.data.frame(rownames(raw))
  expr <- raw
  expr$gene <- rownames(expr)
} else {
  fD <- as.data.frame(sub("\\..*$", "", rownames(raw))) # be sure no trailing number on ensemblId
  expr <- raw
  expr$gene <- sub("\\..*$", "", rownames(expr))
}
#fD        <- as.data.frame(sub("\\..*$", "", rownames(raw))) 
names(fD) <- annoType
head(fD)

# match the gene name with its corresponding biotype
m          <- match(fD[,1], anno[, annoType])
fD$biotype <- anno[m, ]$gene_biotype
head(fD)

# make a new dataframe of gene expression data
# only include counts > 0
# expr      <- raw
# expr$gene <- sub("\\..*$", "", rownames(expr))
expr      <- expr[rowSums(expr[,-ncol(expr)]) > 0, ]

# add geneID and biotype information
m            <- match(expr$gene, fD[,1])
expr$biotype <- fD$biotype[m]
head(expr[, (ncol(expr)-5):ncol(expr)])

# melt for plotting
geneTable <- melt(expr, id = c("gene", "biotype"))
geneTable <- geneTable[geneTable$value > 0, ]

# order biotypes by most --> least expressed
ordered <- as.data.table(geneTable)[, .N, by = biotype][order(-N)]
head(ordered)
#write.csv(ordered, file = paste(dir, "biotypes_ordered.csv", sep = "/"), col.names = TRUE, )

# filter for the first 6 options for cleaner plotting
# all other biotypes are labeled as "other"
filtBio <- as.vector(ordered[6:1,]$biotype)
other   <- as.vector(ordered[7:nrow(ordered)]$biotype)
geneTable$biotype[geneTable$biotype %in% other] <- "other"
geneTable$biotype <- factor(geneTable$biotype, levels = c("other", filtBio))

# sum counts across all biotype/sampleID combinations
sampleTable <- as.data.table(geneTable)[, .N, by = .(variable)]
geneTable   <- as.data.table(geneTable)[, .N, by = .(biotype, variable)]

# get the fraction for each biotype per sample
m              <- match(geneTable$variable, sampleTable$variable)
geneTable$Frac <- geneTable$N / sampleTable[m, ]$N
head(geneTable)

# add in contrast information
#m                  <- match(geneTable$variable, md$SampleID)
m <- match(geneTable$variable, md[[io$sampleID]])
# m <- match(geneTable$variable, md$Alias)
geneTable$contrast <- md[, io$contrast][m]

# keep same order as other plots
levels(geneTable$variable) <- levels(sumCounts.df$SampleID)
head(geneTable)

# barplot of biotypes by sample
biotype.plot <- ggplot(geneTable, aes(x=variable, y=Frac, fill=biotype)) +
  geom_bar(stat="identity") +
  ylab("Fraction of each biotype") +
  xlab("Sample") +
  scale_fill_brewer( palette = "YlGnBu" ) +
  theme(axis.text.x=element_text(angle = 90, size = 8),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
biotype.plot
ggsave(filename = paste(io$outDir, "biotype_barplot.png", sep = "/"), device = "png")


### Plot gene attributes
########################
# create datatable from read distribution file
rDist.dt          <- data.table(rDist, keep.rownames = "sampleID")
rDist.dt$sampleID <- factor(rDist.dt$sampleID, levels = levels(sumCounts.df$SampleID))

# melt for plotting
rDist.melted <- melt(rDist.dt)

# relevel by exon
rDist.melted$variable <- relevel(rDist.melted$variable, ref = "Intergenic")
#rDist.melted$sampleID <- factor(rDist.melted$sampleID, levels = rDist.melted$sampleID)

# barplot of gene fractions by sample
mappings.plot <- ggplot(
  rDist.melted, 
  aes(x = sampleID, y = value, fill = factor(variable, levels = c("Intergenic", "Intron", "Exon")))) +
  geom_bar(position = "fill", stat = "identity") +
  ylab("Fraction of each gene attribute") +
  xlab("Sample") +
  scale_fill_brewer( palette = "Purples" ) +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
mappings.plot
ggsave(filename = paste(io$outDir, "geneAttributes_barplot.png", sep = "/"), device = "png")


### Pairwise sample correlations
################################
# function to save pheatmap
save_pheatmap_png <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, 1000, 2000)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# get metadata columns for heatmap
annoMD <- md
samples <- annoMD[[io$sampleID]]
cols2keep <- unique(c(io$contrast, io$plotCols))

# check dimensions needed for formatting
if (length(cols2keep) >= 2) {
  print("plotting 2+ sample annotations")
  annoMD <- annoMD[, colnames(annoMD) %in% cols2keep]
  annoMD <- annoMD[order(annoMD[, io$contrast], rownames(annoMD)), ]
  ord <- rownames(annoMD)
} else {
  print("plotting 1 sample annotation")
  annoMD <- annoMD[, colnames(annoMD) %in% cols2keep]
  annoMD <- as.data.frame(annoMD, row.names = samples)
  colnames(annoMD) <- cols2keep
}
head(annoMD)

# # factorize annotations
# annoMD$RINcat <- factor(annoMD$RINcat, levels = c("low", "med", "high", NA))

# set the correlation test type
type <- tolower(io$corType)
type

# # define better colors for pheatmap
# lane <- colorRampPalette(c("#1b9e77", "#d95f02", "#7570b3"))(length(unique(annoMD$Lane))) # Dark2
# names(lane) <- unique(annoMD$Lane)[order(unique(annoMD$Lane))]

# group <- colorRampPalette(c("#d7191c", "#ffffbf", "#2c7bb6"))(length(unique(annoMD$Group))) # RdYlBu
# names(group) <- unique(annoMD$Group)[order(unique(annoMD$Group))] 

# ages <- unique(annoMD$Age_at_collection)
# ages <- ages[order(ages)]
# age <- colorRampPalette(c("#ece2f0", "#a6bddb","#1c9099"))(length(ages)) # PuBuGn
# names(age) <- ages

# sex <- colorRampPalette(c("#af8dc3", "#7fbf7b","#f7f7f7"))(length(unique(annoMD$Sex))) # PRGn
# names(sex) <- unique(annoMD$Sex)[order(unique(annoMD$Sex))]

# tm_class <- colorRampPalette(c('#f1a340','#998ec3','#f7f7f7'))(length(unique(annoMD$TM_class))) # PuGr
# names(tm_class) <- unique(annoMD$TM_class)

# novo_rin <- colorRampPalette(c("#9e0142", "#ffffbf", "#5e4fa2"))(length(unique(annoMD$Novogene_RIN))) # Spectral
# names(novo_rin) <- unique(annoMD$Novogene_RIN)[order(unique(annoMD$Novogene_RIN))]

# rin_cat <- colorRampPalette(c("#9e0142", "#ffffbf", "#5e4fa2"))(length(unique(annoMD$RINcat))) # Spectral
# #names(rin_cat) <- unique(annoMD$RINcat)
# names(rin_cat) <- c("low", "med", "high", NA)
# #group_sub <- colorRampPalette(c("#a6cee3", "#1f78b4", "#b2df8a"))(3) # Paired
# #names(group_sub) <- unique(key.dt$patient_dx2)

# # diabetes <- colorRampPalette(c('#f1a340','#998ec3','#f7f7f7'))(3) #PuGr
# # names(diabetes) <- unique(key$diabetes)

# # age_bin <- colorRampPalette(c("#ece2f0", "#1c9099"))(length(unique(annotations$age_bin))) # PuBuGn
# # names(age_bin) <- unique(annotations$age_bin)

# # save colors to list
# ann_colors <- list(
#   Lane = lane,
#   Group = group,
#   Age_at_collection = age,
#   Sex = sex,
#   TM_class = tm_class,
#   Novogene_RIN = novo_rin,
#   RINcat = rin_cat
# )


# for Pearson or Spearman test types
if (type == "pearson" || type == "spearman") {
  print(paste("Running a", type, "correlation test", sep = " "))
  
  # run correlation test and save output
  cor.res <- rcorr(as.matrix(raw), type = type)
  filename <- paste(type, "corrResults.rds", sep = "_")
  saveRDS(cor.res, file = paste(io$outDir, filename, sep= "/"))
  
  # unordered correlation plot
  filename <- paste(type, "unorderedCorrplot_noLabel.png", sep = "_")
  png(filename = paste(io$outDir, filename, sep = "/"), 1000, 1000)
  corrplot(cor.res$r, method = "color", #addCoef.col = "darkgrey", 
           title = paste("Unordered correlation heatmap", io$corType, sep = ": "), 
           outline = TRUE, number.cex = 1.25,
           tl.cex = 1.5, cl.cex = 1.25, mar=c(0,0,2,0))
  dev.off()
  
  # h-clustered correlation plot
  filename <- paste(type, "hclustCorrplot.png", sep = "_")
  png(filename = paste(io$outDir, filename, sep = "/"), 1000, 1000)
  corrplot(cor.res$r, method = "color", order = "hclust",
           title = paste("Hierarchical clustering correlation heatmap", io$corType, sep = ": "),
           number.cex = 1.25, outline = TRUE,
           tl.cex = 0.75, cl.cex = 1.25, mar=c(0,0,2,0))
  dev.off()
  
  # annotated pheatmap
  filename <- paste(type, "annoPheatmap.png", sep = "_")
  plot <- pheatmap(cor.res$r, 
                   main = paste("Clustered sample correlation heatmap", io$corType, sep = ": "),
                   cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D",
                   annotation_col = annoMD, annotation_row = annoMD,
                  #  , annotation_colors = ann_colors,
                   fontsize = 12, fontsize_row = 10, fontsize_col = 10, show_rownames = TRUE, show_colnames = FALSE,
                   treeheight_col = 100,
                   labels_row = as.character(rownames(cor.res$r)),
                   #labels_col = as.character(colnames(cor.res$r)),
                   color = colorRampPalette(c("navy", "white", "red"))(50)
  )
  save_pheatmap_png(plot, paste(io$outDir, filename, sep = "/"))

  # pairwise counts
  sampleNum <- dim(raw)[2]
  samples <- colnames(raw)
  filename <- paste(type, "pairwiseCorrplot.png", sep = "_")
  png(paste(io$outDir, filename, sep = "/"), height = (sampleNum*200), width = (sampleNum*200))
  par(mfrow=c(sampleNum, sampleNum))
  for (i in 1:sampleNum) {
    for (j in 1:sampleNum) {
      plot(log10(raw[,i]), log10(raw[,j]),
           xlab=paste("log10(", paste(samples[i]), ")", sep = ""),
           ylab=paste("log10(", paste(samples[j]), ")", sep = ""),
           col="#1C0DFF15", pch=20, cex.lab = 1.5)
      legend("topleft", legend = round(cor(raw[,i], raw[,j], method = type), 4), cex=1.4)
    }
  }
  dev.off()
}

# run both correlation tests
if (type == "both") {
  print("Running both Pearson and Spearman correlation tests")
  
  ### Pearson correlation
  #######################
  type <- "pearson"
  cor.res <- rcorr(as.matrix(raw), type = type)
  filename <- paste(type, "corrResults.rds", sep = "_")
  saveRDS(cor.res, file = paste(io$outDir, filename, sep= "/"))
  
  # unordered correlation plot
  filename <- paste(type, "unorderedCorrplot.png", sep = "_")
  png(filename = paste(io$outDir, filename, sep = "/"), 1000, 1000)
  corrplot(cor.res$r, method = "color", addCoef.col = "darkgrey", 
           title = paste("Unordered correlation heatmap", "Pearson", sep = ": "), 
           outline = TRUE, number.cex = 1.25,
           tl.cex = 1.5, cl.cex = 1.25, mar=c(0,0,2,0))
  dev.off()
  
  # h-clustered correlation plot
  filename <- paste(type, "hclustCorrplot.png", sep = "_")
  png(filename = paste(io$outDir, filename, sep = "/"), 1000, 1000)
  corrplot(cor.res$r, method = "square", order = "hclust",
           title = paste("Hierarchical clustering correlation heatmap", "Pearson", sep = ": "),
           number.cex = 1.25, outline = TRUE,
           tl.cex = 1.5, cl.cex = 1.25, mar=c(0,0,2,0))
  dev.off()
  
  # annotated pheatmap
  filename <- paste(type, "annoCorrplot.png", sep = "_")
  plot <- pheatmap(cor.res$r, 
           main = paste("Clustered correlation heatmap", "Pearson", sep = ": "),
           cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D",
           annotation_col = annoMD, annotation_row = annoMD, 
           fontsize = 12, fontsize_row = 10, fontsize_col = 10, show_rownames = TRUE,
           labels_row = as.character(rownames(cor.res$r)),
           labels_col = as.character(colnames(cor.res$r)),
           color = colorRampPalette(c("navy", "white", "red"))(50)
  )
  save_pheatmap_png(plot, paste(io$outDir, filename, sep = "/"))
  
  # pairwise counts
  sampleNum <- dim(raw)[2]
  samples <- colnames(raw)
  filename <- paste(type, "pairwiseCorrplot.png", sep = "_")
  png(paste(io$outDir, filename, sep = "/"), height = (sampleNum*200), width = (sampleNum*200))
  par(mfrow=c(sampleNum, sampleNum))
  for (i in 1:sampleNum) {
    for (j in 1:sampleNum) {
      plot(log10(raw[,i]), log10(raw[,j]),
           xlab=paste("log10(", paste(samples[i]), ")", sep = ""),
           ylab=paste("log10(", paste(samples[j]), ")", sep = ""),
           col="#1C0DFF15", pch=20, cex.lab = 1.5)
      legend("topleft", legend = round(cor(raw[,i], raw[,j], method = type), 4), cex=1.4)
    }
  }
  dev.off()
  print("Finished testing with Pearson")
  
  ### Spearman correlation
  ########################
  type <- "spearman"
  cor.res <- rcorr(as.matrix(raw), type = type)
  filename <- paste(type, "corrResults.rds", sep = "_")
  saveRDS(cor.res, file = paste(io$outDir, filename, sep= "/"))

  # unordered correlation plot
  filename <- paste(type, "unorderedCorrplot.png", sep = "_")
  png(filename = paste(io$outDir, filename, sep = "/"), 1000, 1000)
  corrplot(cor.res$r, method = "color", addCoef.col = "darkgrey", 
           title = paste("Unordered correlation heatmap", "Spearman", sep = ": "), 
           outline = TRUE, number.cex = 1.25,
           tl.cex = 1.5, cl.cex = 1.25, mar=c(0,0,2,0))
  dev.off()
  
  # clustered correlation plot
  filename <- paste(type, "hclustCorrplot.png", sep = "_")
  png(filename = paste(io$outDir, filename, sep = "/"), 1000, 1000)
  corrplot(cor.res$r, method = "square", order = "hclust",
           title = paste("Hierarchical clustering correlation heatmap", "Spearman", sep = ": "),
           number.cex = 1.25, outline = TRUE,
           tl.cex = 1.5, cl.cex = 1.25, mar=c(0,0,2,0))
  dev.off()
  
  # annotated pheatmap
  filename <- paste(type, "annoCorrplot.png", sep = "_")
  plot <- pheatmap(cor.res$r, 
           main = paste("Clustered correlation heatmap", "Spearman", sep = ": "),
           cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D",
           annotation_col = annoMD, annotation_row = annoMD,
           fontsize = 12, fontsize_row = 10, fontsize_col = 10, show_rownames = TRUE,
           labels_row = as.character(rownames(cor.res$r)),
           labels_col = as.character(colnames(cor.res$r)),
           color = colorRampPalette(c("navy", "white", "red"))(50)
  )
  save_pheatmap_png(plot, paste(io$outDir, filename, sep = "/"))

  # pairwise counts
  sampleNum <- dim(raw)[2]
  samples <- colnames(raw)
  filename <- paste(type, "pairwiseCorrplot.png", sep = "_")
  png(paste(io$outDir, filename, sep = "/"), 
      height = (sampleNum*200), width = (sampleNum*200))
  par(mfrow=c(sampleNum, sampleNum))
  for (i in 1:sampleNum) {
    for (j in 1:sampleNum) {
      plot(log10(raw[,i]), log10(raw[,j]),
           xlab=paste("log10(", paste(samples[i]), ")", sep = ""),
           ylab=paste("log10(", paste(samples[j]), ")", sep = ""),
           col="#1C0DFF15", pch=20, cex.lab = 1.5)
      legend("topleft", legend = round(cor(raw[,i], raw[,j], method = type), 4), cex=1.4)
    }
  }
  dev.off()
  print("Finished testing with Spearman")
}

# save annoTable to md
# write.table(annoMD, paste(io$outDir, "annoMD.tsv", sep = "/"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
