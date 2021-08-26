# Load libraries
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(scales)

### Set up params
dir <- "/Users/grayaly/CEDAR/PROJECTS"
#annoFile <- "/PLATELET/Best-bulkRNA/BEST2019/CEDAR-pipeline/hg38.Ens_94.biomaRt.geneAnno.Rdata" # Human
#annoFile <- "/ASXL1/ASXL1_seq/CEDAR-pipe-results/mm10.Ens_96.biomaRt.geneAnno.Rdata" # Mouse

# Best 2019 files
# metaFile <- "/PLATELET/Best-bulkRNA/BEST2019/DATA/metadata.tsv"
# countsFile <- "/PLATELET/Best-bulkRNA/BEST2019/CEDAR-pipeline/Best_GSE107868_counts.txt"
# readDistFile <- "PLATELET/Best-bulkRNA/BEST2019/CEDAR-pipeline/results/tables/read_coverage.txt"

# Best 2015 files
# metaFile <- "/PLATELET/Best-bulkRNA/BEST2015/CEDAR-pipeline/metadata.tsv"
# countsFile <- "/PLATELET/Best-bulkRNA/BEST2015/CEDAR-pipeline/Best_GSE68086_counts.txt"
# readDistFile <- "PLATELET/Best-bulkRNA/BEST2015/CEDAR-pipeline/results/tables/read_coverage.txt"

# MaxsonLab ASXL1 files
# metaFile <- "/ASXL1/ASXL1_seq/CEDAR-pipe-results/metadata.tsv"
# countsFile <- "/ASXL1/ASXL1_seq/CEDAR-pipe-results/MaxsonLab_ASXL1_counts.txt"
# rDistFile <- "/ASXL1/ASXL1_seq/CEDAR-pipe-results/tables/read_coverage.txt"

### Load and organize data
# Read in annotation data
anno <- get(load(paste(dir, annoFile, sep = "/")))

# Read in project data; organize sampleIDs so each are the same
md <- read.table(paste(dir, metaFile, sep = "/"), stringsAsFactors = FALSE, sep = "\t", header = TRUE)
#colnames(md) <- tolower(colnames(md))
#md <- md[order(md$Group, md$Rep), ] ## Need to make agnostic
ord <- md$SampleID #### Check capital!!!

raw <- read.table(paste(dir, countsFile, sep = "/"), header = TRUE, sep = "\t", row.names = 1)
raw <- raw[, colnames(raw) %in% md$SampleID]
raw <- raw[, ord]

rDist <- read.table(paste(dir, rDistFile, sep = "/"), header = TRUE, sep = "\t", row.names = 1)
rDist <- rDist[rownames(rDist) %in% md$SampleID, ]
rDist <- rDist[ord, ]

### Plot sumCounts
sumCounts.df <- as.data.frame(apply(raw, 2, sum))
names(sumCounts.df) <- "sumCounts"
sumCounts.df$SampleID <- rownames(sumCounts.df)
sumCounts.df$SampleID <- factor(sumCounts.df$SampleID, levels = sumCounts.df$SampleID) # keeps the correct numerical order for ggplot
sumCounts.df$contrast <- md$Condition ### Write differently for each metadata comparison

# Barplot of read sumCounts for each sample
totalReads.plot <- ggplot(
  sumCounts.df, 
  aes(x = SampleID, y = sumCounts, fill = contrast)) +
  geom_col() +
  #ggtitle("Total Reads per Sample") +
  ylab("Number of Input Reads") +
  xlab("Sample") +
  scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(labels = scientific) +
  #scale_fill_brewer(palette = "YlGnBu") +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
#totalReads.plot

### Plot biotypes
# Check how gene names are annotated
annoType <- c()
if (is.null(grep("ENSG", rownames(raw)))) {
  annoType <- "ensembl_gene_id"
} else {
  annoType <- "external_gene_name"
}

# Get gene names from counts matrix
fD <- as.data.frame(rownames(raw))
names(fD) <- annoType

# Match the gene name with its corresponding biotype
m <- match(fD[,1], anno[, annoType])
fD$biotype <- anno[m, ]$gene_biotype

# Make new DF of gene expression data; only include counts > 0
expr <- raw
expr$gene <- rownames(expr)
expr <- expr[rowSums(expr[,-ncol(expr)]) > 0, ]

# Add geneID and biotype information
m <- match(expr$gene, fD[,1])
expr$biotype <- fD$biotype[m]

# Melt for plotting
geneTable <- melt(expr, id = c("gene", "biotype"))
geneTable <- geneTable[geneTable$value > 0, ]

# Order biotypes by most --> least expressed
ordered <- as.data.table(geneTable)[, .N, by = biotype][order(-N)]
#write.csv(ordered, file = paste(dir, "biotypes_ordered.csv", sep = "/"), col.names = TRUE, )

# Filter for the first 6 options for cleaner plotting; all other biotypes are labeled as "other"
filtBio <- as.vector(ordered[6:1,]$biotype)
other <- as.vector(ordered[7:nrow(ordered)]$biotype)
geneTable$biotype[geneTable$biotype %in% other] <- "other"
geneTable$biotype <- factor(geneTable$biotype, levels = c("other", filtBio))

# Sum counts across all biotype/sampleID combinations
sampleTable <- as.data.table(geneTable)[, .N, by = .(variable)]
geneTable <- as.data.table(geneTable)[, .N, by = .(biotype, variable)]

# Get the fraction for each biotype per sample
m <- match(geneTable$variable, sampleTable$variable)
geneTable$Frac <- geneTable$N / sampleTable[m, ]$N

# Add in contrast information
m <- match(geneTable$variable, md$SampleID)
geneTable$contrast <- md$Condition[m] ### Change to snakemake[contrast] ###

# Plot biotype bar
biotype.plot <- ggplot(geneTable, aes(x=variable, y=Frac, fill=biotype)) +
  geom_bar(stat="identity") +
  ylab("Fraction of each biotype") +
  xlab("Sample") +
  scale_fill_brewer( palette = "YlGnBu" ) +
  theme(axis.text.x=element_text(angle = 90, size = 5),
        axis.text.y=element_text(size=5),
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
#biotype.plot

### Plot gene attributes
rDist.dt <- data.table(rDist, keep.rownames = "sampleID")
rDist.dt$sampleID <- factor(rDist.dt$sampleID, levels = rDist.dt$sampleID)

# Stacked barplot
rDist.melted <- melt(rDist.dt)

# relevel by exon
rDist.melted$variable <- relevel(rDist.melted$variable, ref = "Intergenic")
#rDist.melted$sampleID <- factor(rDist.melted$sampleID, levels = rDist.melted$sampleID)

mappings.plot <- ggplot(
  rDist.melted, 
  aes(x = sampleID, y = value, fill = factor(variable, levels = c("Intergenic", "Intron", "Exon")))) +
  geom_bar(position = "fill", stat = "identity") +
  ylab("Fraction of each gene attribute") +
  xlab("Sample") +
  scale_fill_brewer( palette = "Purples" ) +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
#mappings.plot

### Arrange plots together
arranged.plot <- ggarrange(totalReads.plot, mappings.plot, biotype.plot,
                           ncol = 1, nrow = 3, labels = "AUTO", align = "h")
#print(arranged.plot)
