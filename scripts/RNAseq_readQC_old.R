
### Set up commandArgs
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
  cat("--outdir       : Path to the directory to save the figures in.           [ required ]
      \n")
  cat("\n")
  q()
}
## Save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  annoFile     <- sub('--annoFile=', '', args[grep('--annoFile=', args)])
  metaFile     <- sub('--metaFile=', '', args[grep('--metaFile=', args)])
  countsFile   <- sub('--countsFile=', '', args[grep('--countsFile=', args)])
  readDistFile <- sub('--readDistFile=', '', args[grep('--readDistFile=', args)])
  contrast     <- sub('--contrast=', '', args[grep('--contrast=', args)])
  outDir       <- sub('--outdir=', '', args[grep('--outdir=', args)])
}

# create outdir as needed
if(!(file.exists( outDir ))) {
    print(paste("mkdir:", outDir))
    dir.create(outDir, FALSE, TRUE)  
}

## Check input files in logs/.out file
io <- list(
  annoFile = annoFile,
  metaFile = metaFile,
  countsFile = countsFile,
  readDistFile = readDistFile,
  contrast = contrast,
  outDir = outDir
)
io

# Debugging
#io <- list(
#  annoFile = "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_94.biomaRt.geneAnno.Rdata",
#  metaFile = "data/PP_metadata_keep_FINAL_updated.txt",
#  countsFile = "data/PP_cohort_Run3_counts.txt",
#  readDistFile = "results/tables/read_coverage.txt",
#  contrast = "Status",
#  outDir = "results/readQC_plots"
#)
#io

### Load libraries
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(scales)

### Load and organize data
# Read in annotation data
load(io$annoFile, verbose = TRUE)

# Read in meta data; organize sampleIDs so the order is the same
md  <- read.table(io$metaFile, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
md  <- md[order(md[, io$contrast], md[, "SampleID"]) , ]
ord <- md$SampleID

# Read in and reorder unfiltered counts table
raw <- read.table(io$countsFile, sep = "\t", header = TRUE, row.names = 1)
raw <- raw[, colnames(raw) %in% md$SampleID]
raw <- raw[, ord]

# Read in and reorder read_coverage.txt
rDist <- read.table(io$readDistFile, header = TRUE, sep = "\t", row.names = 1)
rDist <- rDist[rownames(rDist) %in% md$SampleID, ]
rDist <- rDist[ord, ]

### Plot sumCounts
sumCounts.df          <- as.data.frame(apply(raw, 2, sum))
names(sumCounts.df)   <- "sumCounts"
sumCounts.df$SampleID <- rownames(sumCounts.df)
iv                    <- match(rownames(sumCounts.df), md$SampleID)
sumCounts.df$contrast <- md[iv, io$contrast]
stopifnot(rownames(sumCounts.df) == md[iv,"SampleID"])

sumCounts.df$SampleID <- factor(sumCounts.df$SampleID, 
                                levels = paste(sumCounts.df[order(sumCounts.df$contrast, -sumCounts.df$sumCounts),"SampleID"])) # keeps the correct numerical order for ggplot
sumCounts.df$contrast <- md[, io$contrast]
ord <- levels(sumCounts.df$SampleID)

# Barplot of read sumCounts for each sample
totalReads.plot <- ggplot(
    sumCounts.df, 
    aes(x = SampleID, y = sumCounts, fill = contrast)) +
    geom_col(color="black") +
    ylab("Number of input reads") +
    xlab("Sample") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scientific) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = "black"))

# Violin plot of read sumCounts
totalReads.boxplot <- ggplot(
    sumCounts.df,
    aes(x = contrast, y = sumCounts, fill = contrast)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(color = "black", width = 0.1, fill = "white") +
    xlab("Contrast") +
    ylab("Numer of input reads") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(labels = scientific) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(colour = "black"))

### Plot mitochondrial fraction
# Check how gene names are annotated in counts
annoType <- c()
if( length(grep("ENSG", rownames(raw))) == 0) {
  print("counts are external gene ids")
  annoType <- "external_gene_name"
  mtSub.mat             <- raw[grep("MT-", rownames(raw), value = TRUE, ignore.case = TRUE), ]
  sumCounts.df$mtCounts <- colSums(mtSub.mat)
  sumCounts.df$fracMT   <- sumCounts.df$mtCounts / sumCounts.df$sumCounts
}else{
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

# Barplot of mito read fraction
mtReads.plot <- ggplot(
    sumCounts.df, 
    aes(x = SampleID, y = fracMT, fill = contrast)) +
    geom_col(color="black") +
    ylab("Percent mitochondrial reads") +
    xlab("Sample") +
    scale_fill_brewer(palette = "Paired") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(colour = "black"))

# Violinplot of mito read fraction
mtReads.boxplot <- ggplot(
    sumCounts.df,
    aes(x = contrast, y = fracMT, fill = contrast)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(color = "black", width = 0.1, fill = "white") +
    xlab("Contrast") +
    ylab("Percent mitochondrial reads") +
    scale_fill_brewer(palette = "Paired") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.line = element_line(colour = "black"))


# Save summary plot
pdf(paste(io$outDir, "readSummary.pdf", sep = "/"), onefile = FALSE)
ggarrange(totalReads.plot, totalReads.boxplot,
          mtReads.plot, mtReads.boxplot,
          ncol = 2, nrow = 2, labels = "AUTO", common.legend = TRUE, legend = "right")
dev.off()

png(paste(io$outDir, "readSummary.png", sep = "/"), 1000, 1000)
ggarrange(totalReads.plot, totalReads.boxplot,
          mtReads.plot, mtReads.boxplot,
          ncol = 2, nrow = 2, labels = "AUTO", common.legend = TRUE)
dev.off()

### Plot biotypes

# Get gene names from counts matrix
if( !annoType == "ensembl_gene_id") {
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

# Match the gene name with its corresponding biotype
m          <- match(fD[,1], anno[, annoType])
fD$biotype <- anno[m, ]$gene_biotype

# Make new DF of gene expression data; only include counts > 0
# expr      <- raw
# expr$gene <- sub("\\..*$", "", rownames(expr))
expr      <- expr[rowSums(expr[,-ncol(expr)]) > 0, ]

# Add geneID and biotype information
m            <- match(expr$gene, fD[,1])
expr$biotype <- fD$biotype[m]

# Melt for plotting
geneTable <- melt(expr, id = c("gene", "biotype"))
geneTable <- geneTable[geneTable$value > 0, ]

# Order biotypes by most --> least expressed
ordered <- as.data.table(geneTable)[, .N, by = biotype][order(-N)]
#write.csv(ordered, file = paste(dir, "biotypes_ordered.csv", sep = "/"), col.names = TRUE, )

# Filter for the first 6 options for cleaner plotting; all other biotypes are labeled as "other"
filtBio <- as.vector(ordered[6:1,]$biotype)
other   <- as.vector(ordered[7:nrow(ordered)]$biotype)

geneTable$biotype[geneTable$biotype %in% other] <- "other"
geneTable$biotype <- factor(geneTable$biotype, levels = c("other", filtBio))

# Sum counts across all biotype/sampleID combinations
sampleTable <- as.data.table(geneTable)[, .N, by = .(variable)]
geneTable   <- as.data.table(geneTable)[, .N, by = .(biotype, variable)]

# Get the fraction for each biotype per sample
m              <- match(geneTable$variable, sampleTable$variable)
geneTable$Frac <- geneTable$N / sampleTable[m, ]$N

# Add in contrast information
m                  <- match(geneTable$variable, md$SampleID)
geneTable$contrast <- md[, io$contrast][m]

# keep same order as other plots
levels(geneTable$variable) <- levels(sumCounts.df$SampleID)

# Plot biotype bar
pdf(paste(io$outDir, "biotypes.pdf", sep = "/"), width=7, height=5)
biotype.plot <- ggplot(geneTable, aes(x=variable, y=Frac, fill=biotype)) +
  geom_bar(stat="identity") +
  ylab("Fraction of each biotype") +
  xlab("Sample") +
  scale_fill_brewer( palette = "YlGnBu" ) +
  theme(axis.text.x=element_text(angle = 90, size = 5),
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
dev.off()

### Plot gene attributes
rDist.dt          <- data.table(rDist, keep.rownames = "sampleID")
rDist.dt$sampleID <- factor(rDist.dt$sampleID, levels = levels(sumCounts.df$SampleID))

# Melt for plotting
rDist.melted <- melt(rDist.dt)
# relevel by exon
rDist.melted$variable <- relevel(rDist.melted$variable, ref = "Intergenic")
#rDist.melted$sampleID <- factor(rDist.melted$sampleID, levels = rDist.melted$sampleID)


pdf(paste(io$outDir, "mappings.pdf", sep = "/"),width=7, height=5)
mappings.plot <- ggplot(
  rDist.melted, 
  aes(x = sampleID, y = value, fill = factor(variable, levels = c("Intergenic", "Intron", "Exon")))) +
  geom_bar(position = "fill", stat = "identity") +
  ylab("Fraction of each gene attribute") +
  xlab("Sample") +
  scale_fill_brewer( palette = "Purples" ) +
  theme(axis.text.x = element_text(angle = 90, size = 5),
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
dev.off()

