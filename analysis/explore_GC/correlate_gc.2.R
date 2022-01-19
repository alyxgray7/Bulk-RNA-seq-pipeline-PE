########################################################
#####   Script to correlate %GC values to other    #####
#####           known sample variables.            #####
########################################################

# Summarizes %GC distributions in supplemental plots.
# Correlates %GC distributions with other known sample variables to look for patterns.

# Inputs: Compiled %GC distributions by sample.
# Outputs:


########################################################
#####                   Set up                     #####
########################################################

### Libraries
#############
library(ggplot2)
library(data.table)
library(scales)
library(RColorBrewer)
# library(Hmisc)
# library(pheatmap)
# library(corrplot)


### Std in/out
##############
io <- list(
    # gcFile = "results/tables/compiled_readGC_test.tsv",
    gcFile = "results/tables/compiled_readGC.tsv",
    # mdFile = "data/metadata.tsv",
    mdFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/results/starSummary/md_merged.tsv",
    outDir = "analysis/explore_GC",
    sampleID = "rna_sampleID"
)
io

# create outdir as needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}


### Functions
#############
# function to save pheatmap
save_pheatmap_png <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, 1000, 2000)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


### Format data
###############
# metadata
md <- read.table(io$mdFile, header = TRUE, row.names = NULL, sep = "\t")
md$rna_sampleID <- tstrsplit(md$sampleID, split = "_", fixed = TRUE)[[1]]
head(md[[io$sampleID]])

# gc distribution file
gc <- read.table(io$gcFile, header = TRUE, row.names = 1, sep = "\t")
colnames(gc) <- tstrsplit(colnames(gc), split = "_", fixed = TRUE)[[1]]
gc$pct_GC <- rownames(gc)
head(gc)
dim(gc)


########################################################
#####           Summarize GC Distributions         #####
########################################################

### Summarize GC
################
# melt table for ggplot
gc.melted <- reshape2::melt(gc, value.name = "count")
gc.melted$variable <- as.factor(gc.melted$variable)
gc.melted$pct_GC <- as.numeric(gc.melted$pct_GC)
gc.melted <- na.omit(gc.melted)

# add other variables to melted table
m <- match(gc.melted$variable, md$rna_sampleID)
gc.melted$RINcat <- md$RINcat[m]
gc.melted$RINcat <- factor(gc.melted$RINcat, levels = c("low", "med",'high'))
gc.melted$RIN <- md$Novogene_RIN[m]
gc.melted$pctUnmapped <- md$pctUnmapped[m]
gc.melted$numMappedReads <- md$numUniqMappedReads[m]
head(gc.melted)

# get the mean pct_GC for each sample
gcsMean <- list()
for (mysample in levels(gc.melted$variable)) {
    sub <- gc.melted[which(gc.melted$variable == mysample), ]
    sub$prod <- sub$pct_GC * sub$count
    gcsMean[[mysample]] <- ( sum(sub$prod) / sum(sub$count) )
}
gcsMean.melted <- reshape2::melt(as.data.frame(gcsMean))
colnames(gcsMean.melted) <- c("sampleID", "mean_pctGC")
head(gcsMean.melted)

# get the highest frequency pct_GC for each sample
gcPeak <- list()
gcsFreq <- apply(gc[, -which(colnames(gc) == "pct_GC")], 2, function(x) max(x, na.rm = TRUE))
for (myname in names(gcsFreq)) {
    gcPeak[[myname]] <- gc[["pct_GC"]][which(gc[[myname]] == max(gc[[myname]], na.rm = TRUE))]
}
gcPeak.mat <- as.matrix(gcPeak)

# merge with metadata
md <- merge(md, gcsMean.melted, by.x = sampleID, by.y = "sampleID")
md$gcPeak <- rep(NA, nrow(md))
for (i in 1:nrow(gcPeak.mat)) {
    md[i, "gcPeak"] <- gcPeak.mat[,1][[i]]
}
head(md)
write.table(md, paste(io$outDir, "merged_md.tsv", sep = "/"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


### Summary plots
#################
# facet histogram, colored by sample
ggplot(gc.melted, aes(x = pct_GC, y = count, color = variable)) +
    geom_col(show.legend = FALSE) +
    geom_vline(aes(xintercept = value), gcsMean.melted, color = "red", linetype = "dashed") +
    facet_wrap(~ variable) +
    xlab("Percent GC") +
    ylab("Count") +
    theme_bw() +
    #theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, size = 8), 
          axis.text.y = element_text(size = 8),
          legend.position = "None")
ggsave(paste(io$outDir, "gcDist_facet.png", sep = "/"), device = "png")

# density plot colored by sample
ggplot(gc.melted, aes(x = pct_GC, color = variable)) +
    geom_density() +
    xlab("Percent GC") +
    ylab("Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 8), 
          axis.text.y = element_text(size = 8),
          legend.position = "None")
ggsave(paste(io$outDir, "gcDensity_bySample.png", sep = "/"), device = "png")

# density plots colored by RINcat
ggplot(gc.melted, aes(x = pct_GC, color = RINcat)) +
    geom_density() +
    xlab("Percent GC") +
    ylab("Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 8), 
          axis.text.y = element_text(size = 8),
          legend.position = "right")
ggsave(paste(io$outDir, "gcDensity_byRINcat.png", sep = "/"), device = "png")

# facet grid by sample
ggplot(gc.melted, aes(x = pct_GC, color = variable)) +
    geom_density() +
    facet_wrap(~variable) +
    geom_vline(aes(xintercept = value), gcsMean.melted, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "None")
ggsave(paste(io$outDir, "gcDensity_facet.png", sep = "/"), device = "png")

# barplot of gcPeaks, colored by RINcat
ggplot(md, aes_string(x = sampleID, y = "gcPeak", fill = "RINcat")) +
    geom_col() +
    theme_bw() +
    xlab("RNA SampleID") +
    ylab("Peak %GC") +
    theme(axis.text.x = element_text(angle = 90, size = 6))
ggsave(paste(io$outDir, "gcPeak_bySample_RINcat.png", sep = "/"), device = "png")

# barplot of gcPeaks, ordered by pctUnmapped
toplot <- md
toplot <- toplot[order(toplot$pctUnmapped), ]
toplot[[sampleID]] <- factor(toplot[[sampleID]], levels = toplot[[sampleID]])
ggplot(toplot, aes_string(x = sampleID, y = "gcPeak")) +
    geom_col() +
    theme_bw() +
    xlab("RNA SampleID (in decreasing pctUnmapped)") +
    ylab("Peak %GC") +
    theme(axis.text.x = element_text(angle = 90, size = 6))
ggsave(paste(io$outDir, "gcPeak_bySample_orderedPctUnmapped.png", sep = "/"), device = "png")


########################################################
#####             Correlate Mean %GC               #####
########################################################

### Pearson correlation test
############################
# format table for correlation test
cols2keep <- c(sampleID, "avgMappedReadLength", "Novogene_RIN", "pctUnmapped", "mean_pctGC", "gcPeak")
test <- md[, which(colnames(md) %in% cols2keep)]
rownames(test) <- test[[sampleID]]
test <- test[, -which(colnames(test) == sampleID)]
head(test)

# apply Pearson correlation test to other variables
cor.test <- rcorr(as.matrix(test), type = "pearson")


### Correlation plots
#####################
# heatmap of all coefficients
png(paste(io$outDir, "heatmap_pearson.png", sep= "/"), 1000, 1000)
corrplot(cor.test$r, method = "color", addCoef.col = "darkgrey", 
           #p.mat = cor.test$P,
           title = paste("Pearson correlation heatmap"), 
           outline = TRUE, number.cex = 1.25,
           tl.cex = 1.5, cl.cex = 1.25, mar=c(0,0,2,0))
dev.off()

# heatmap of significant p-vals
png(paste(io$outDir, "heatmap_pearsonPval.png", sep= "/"), 1000, 1000)
corrplot(cor.test$r, p.mat = cor.test$P,
        method = "color", 
        type = "lower", insig = "blank",
        addCoef.col = "darkgrey", 
        diag = FALSE,
        title = paste("Pearson correlation heatmap (pval < 0.05)"), 
        outline = TRUE, number.cex = 1.25,
        tl.cex = 1.5, cl.cex = 1.25, mar=c(0,0,2,0))
dev.off()

# dotplots of significant relationships
# function to plot
# plotDotplot <- function(data, x, y) {
#     formula <- as.formula(paste(x, y, sep = "~"))
#     r2 <- summary(lm(formula, data = data))$r.squared
#     p <- ggplot(data, aes_string(x = x, y = y)) +
#         geom_point() +
#         geom_smooth(method = "lm")
#     return(print(p))
#     # print(p)
#     # ggsave(paste(io$outDir, filename, sep = "/"), device = "png")
# }

# avgMappedReadLength ~ pctUnmapped
x <- "pctUnmapped"
y <- "avgMappedReadLength"
filename <- paste(paste(x, y, sep = "_"), "png", sep = ".")
r2 <- summary(lm(as.formula(paste(y, x, sep = "~")), data = md))$r.squared
ggplot(md, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    annotate(geom = "text", x = 75, y = 250, label = paste0("R2 = ", r2), color = "red")
ggsave(paste(io$outDir, filename, sep = "/"), device = "png")

# avgMappedReadLength ~ mean_pctGC
x <- "mean_pctGC"
y <- "avgMappedReadLength"
filename <- paste(paste(x, y, sep = "_"), "png", sep = ".")
r2 <- summary(lm(as.formula(paste(y, x, sep = "~")), data = md))$r.squared
ggplot(md, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    annotate(geom = "text", x = 50, y = 250, label = paste0("R2 = ", r2), color = "red")
ggsave(paste(io$outDir, filename, sep = "/"), device = "png")

# avgMappedReadLength ~ gcPeak
x <- "gcPeak"
y <- "avgMappedReadLength"
filename <- paste(paste(x, y, sep = "_"), "png", sep = ".")
r2 <- summary(lm(as.formula(paste(y, x, sep = "~")), data = md))$r.squared
ggplot(md, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    annotate(geom = "text", x = 5, y = 250, label = paste0("R2 = ", r2), color = "red")
ggsave(paste(io$outDir, filename, sep = "/"), device = "png")

# Novogene_RIN ~ pctUnmapped
x <- "pctUnmapped"
y <- "Novogene_RIN"
filename <- paste(paste(x, y, sep = "_"), "png", sep = ".")
r2 <- summary(lm(as.formula(paste(y, x, sep = "~")), data = md))$r.squared
ggplot(md, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    annotate(geom = "text", x = 50, y = 8, label = paste0("R2 = ", r2), color = "red")
ggsave(paste(io$outDir, filename, sep = "/"), device = "png")

# pctUnmapped ~ mean_pctGC
x <- "pctUnmapped"
y <- "mean_pctGC"
filename <- paste(paste(x, y, sep = "_"), "png", sep = ".")
r2 <- summary(lm(as.formula(paste(y, x, sep = "~")), data = md))$r.squared
ggplot(md, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm") +
    annotate(geom = "text", x = 75, y = 45, label = paste0("R2 = ", r2), color = "red")
ggsave(paste(io$outDir, filename, sep = "/"), device = "png")

# pctUnmapped ~ gcPeak
x <- "pctUnmapped"
y <- "gcPeak"
filename <- paste(paste(x, y, sep = "_"), "png", sep = ".")
r2 <- summary(lm(as.formula(paste(y, x, sep = "~")), data = md))$r.squared
ggplot(md, aes_string(x = x, y = y)) +
    geom_point() +
    #geom_smooth(method = "lm") 
    annotate(geom = "text", x = 75, y = 45, label = paste0("R2 = ", r2), color = "red")
ggsave(paste(io$outDir, filename, sep = "/"), device = "png")

# mean_pctGC ~ gcPeak
x <- "mean_pctGC"
y <- "gcPeak"
filename <- paste(paste(x, y, sep = "_"), "png", sep = ".")
r2 <- summary(lm(as.formula(paste(y, x, sep = "~")), data = md))$r.squared
ggplot(md, aes_string(x = x, y = y)) +
    geom_point() +
    #geom_smooth(method = "lm") 
    annotate(geom = "text", x = 75, y = 45, label = paste0("R2 = ", r2), color = "red")
ggsave(paste(io$outDir, filename, sep = "/"), device = "png")
