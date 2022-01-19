########################################################################
#####         Script to summarize STAR mapping results by          #####
####                     metadata information.                     #####
########################################################################

### Set up
##########
# command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
help <- function(){
  cat("sumMapping.R :
  
  - Outputs :
  
  1) 
      ")
  cat("\n")
  cat("Usage : \n")
  cat("--starSummary  : Compiled STAR statistics                                [ required ]
                           (results/tables/{project_id}_STAR_mapping_statistics.txt
      \n")
  cat("--metaFile     : Path to config['omic_meta_data']                        [ required ]
      \n")
  cat("--plotCols     : Column names in metadata file to annotate clustering;   [ required ]
                        specified by config['meta_columns_to_plot']
      \n")
  cat("--outdir       : Path to the directory to save the figures in.           [ required ]
      \n")
  cat("\n")
  q()
}

# save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  starSummary     <- sub('--starSummary=', '', args[grep('--starSummary=', args)])
  metaFile     <- sub('--metaFile=', '', args[grep('--metaFile=', args)])
  plotCols     <- sub('--plotCols=', '', args[grep('--plotCols=', args)])
  outDir       <- sub('--outdir=', '', args[grep('--outdir=', args)])
}

# make plotCols argument R-readable
plotCols <- c(strsplit(plotCols, split = ",", fixed = TRUE)[[1]])

# check input files in logs/.out file
io <- list(
  starSummary = starSummary,
  metaFile = metaFile,
  plotCols = plotCols,
  outDir = outDir
)
io

# debugging on exa
io <- list(
  starSummary = "results/tables/platelet_full-cohort_STAR_mapping_statistics.txt",
  metaFile = "data/metadata.tsv",
  plotCols = c("Lane", "Group", "Sex", "Age_at_collection", "Novogene_RIN", "TM_class"),
  outDir = "results/starSummary"
)
io

# create outdir as needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}

# libraries
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(reshape2)
library(scales)


### Load data
#############
# read in metadata table
md <- as.data.table(read.table(io$metaFile, sep = "\t", header = TRUE, row.names = 1),keep.rownames = "Alias")

# split RNA sample identifiers
md[, "rnaSampleID" := tstrsplit(Alias, split = "_", fixed = TRUE)[[1]]]
md[, "rnaBatch" := gsub("*[0-9]", "", rnaSampleID)]
md[, "rnaBatchOrder" := gsub("*[A-Z]", "", rnaSampleID)]

# bin RIN values
md[Novogene_RIN <= 10.0, "RINcat" := "high"]
md[Novogene_RIN <= 6.0, "RINcat" := "med"]
md[Novogene_RIN <= 4.0, "RINcat" := "low"]

# bin age groups
md[, "Age_cat" := cut(Age_at_collection, 3)]

# factorize
md$rnaBatch <- as.factor(md$rnaBatch)
md$RINcat <- factor(md$RINcat, levels = c("low", "med", "high"))

head(md)
dim(md)


### STAR summary report
#######################
# read in table
star <- read.table(io$starSummary, sep = "\t", header = TRUE, row.names = 1)
head(star[, 1:5])
dim(star)

# edit column names
colnames(star) <- gsub("_bam", "", colnames(star))

# select only samples included in metadata file
star <- star[, which(colnames(star) %in% md$Alias) ]
dim(star)


### Create mapping results table
################################
# extract relevant results
star_sub <- star[c(
  grep("Number of input reads", rownames(star)),
  grep("Average input read length", rownames(star)),
  grep("Uniquely mapped reads number", rownames(star)),
  grep("Uniquely mapped reads %", rownames(star)),
  grep("Average mapped length", rownames(star)),
  grep("% of reads unmapped: too many mismatches", rownames(star)),
  grep("% of reads unmapped: too short", rownames(star)),
  grep("% of reads unmapped: other", rownames(star))
), ]
dim(star_sub)

# clean up table
rownames(star_sub) <- c(
  "numInputReads", 
  "avgInputReadLength", 
  "numUniqMappedReads", 
  "pctUniqMappedReads", 
  "avgMappedReadLength", 
  "pctReadsUnmapped.mismatch", 
  "pctReadsUnmapped.tooshort", 
  "pctReadsUnmapped.other"
)
star_sub <- as.data.frame(t(as.matrix(star_sub)))
ord <- rownames(star_sub)
star_sub <- apply(star_sub, 2, function(x) gsub("%", "", x))
star_sub <- apply(star_sub, 2, as.numeric)
star_sub <- as.data.frame(star_sub)
star_sub$sampleID <- ord
head(star_sub)

# add metadata info
merged <- merge(star_sub, as.data.frame(md), by.x = "sampleID", by.y = "Alias")
merged$pctUnmapped <- (
  merged$pctReadsUnmapped.tooshort +
  merged$pctReadsUnmapped.mismatch +
  merged$pctReadsUnmapped.other
)
write.table(merged, paste(io$outDir, "md_merged.tsv", sep = "/"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


### Summary plots
#################
# color and order by plotCols
plotCols <- io$plotCols
plotCols[grep("Novogene_RIN", plotCols)] <- "RINcat"
plotCols[grep("Age_at_collection", plotCols)] <- "Age_cat"
plotCols <- append(plotCols, "rnaBatch", length(plotCols))
cols2plot <- c("avgMappedReadLength", "pctUniqMappedReads", "pctUnmapped")
for (myfill in plotCols) {
  # myfill <- plotCols[length(plotCols)]
  for (mycol in cols2plot) {
    # mycol <- cols2plot[length(cols2plot)]

    # reorder samples for a cleaner plot
    toplot <- merged
    toplot <- toplot[order(toplot[, myfill], toplot[, mycol]), ]
    toplot$sampleID <- factor(toplot$sampleID, levels = toplot$sampleID)
    if (myfill == "RINcat") {
      toplot[, myfill] <- factor(toplot[, fill], levels = c("low", "med", "high"))
    } else {
      toplot[, myfill] <- as.factor(toplot[, myfill])
    }

    # set palette for plots
    colorCount <- length(unique(toplot[, myfill]))
    getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

    # barplot by sample
    p <- ggplot(toplot, aes_string(x = "sampleID", y = mycol, fill = myfill)) +
      geom_col() +
      geom_col(color = "black") +
      xlab("Sample") +
      #scale_fill_brewer(palette = "Paired") +
      scale_fill_manual(values = getPalette(colorCount)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        # legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = "black"))
    print(p)
    ggsave(paste(io$outDir, paste(myfill, mycol, "barplot.png", sep = "_"), sep = "/"), device = "png")

    # violin + box plot by plotCol
    cols2keep <- which(colnames(toplot) %in% c("sampleID", myfill, mycol))
    toplot <- toplot[, c(cols2keep)]
    toplot.melted <- reshape2::melt(toplot, id.vars = c("sampleID", myfill), value.name = mycol)
    p <- ggplot(toplot, aes_string(x = myfill, y = mycol, fill = myfill)) +
          geom_violin(trim = FALSE) +
          geom_boxplot(color = "black", width = 0.1, fill = "white") +
          # scale_fill_brewer(palette = "Paired") +
          scale_fill_manual(values = getPalette(colorCount)) +
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, size = 8),
                axis.text.y = element_text(size = 10),
                axis.title.y = element_text(size = 10),
                axis.title.x = element_text(size = 10),
                # legend.title = element_blank(),
                legend.text = element_text(size = 10),
                axis.line = element_line(colour = "black"))
    print(p)
    ggsave(paste(io$outDir, paste(myfill, mycol, "boxplot.png", sep = "_"), sep = "/"), device = "png")
  }
}

# boxplot of pctUnmapped by RINcat
ggplot(merged, aes(x = RINcat, y = pctUnmapped)) +
  geom_boxplot()
ggsave(paste(io$outDir, "RINcat_pctUnmapped.png", sep = "/"), device = "png")

# boxplot and dotplot by rnaBatch
ggplot(merged, aes(x = rnaBatch, y = pctUnmapped, color = RINcat)) +
  geom_boxplot(color = "black") +
  geom_point() +
  scale_color_brewer(palette = "Set1")
ggsave(paste(io$outDir, "rnaBatch_pctUnmapped.png", sep = "/"), device = "png")

# barplot of pctUnmapped by sample
# x axis is ordered by RINcat
# colored by RNAbatch
toplot <- merged
toplot <- toplot[order(toplot$Novogene_RIN), ]
toplot$sampleID <- factor(toplot$sampleID, levels = toplot$sampleID)
head(toplot)
colorCount <- length(levels(toplot$rnaBatch))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
ggplot(toplot, aes(x = sampleID, y = pctUnmapped, fill = rnaBatch)) +
  geom_col() +
  scale_fill_manual(values = getPalette(colorCount)) +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        # legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = "black"))
ggsave(paste(io$outDir, "test3.png", sep = "/"), device = "png")

# colored by RINcat
colorCount <- length(levels(toplot$RINcat))
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
ggplot(toplot, aes(x = sampleID, y = pctUnmapped, fill = RINcat)) +
  geom_col() +
  scale_fill_manual(values = getPalette(colorCount)) +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        # legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = "black"))
ggsave(paste(io$outDir, "test4.png", sep = "/"), device = "png")

# density plot of pctUnmapped

