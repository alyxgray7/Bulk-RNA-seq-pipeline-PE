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
library(Hmisc)
library(languageserver)


### Std in/out
##############
io <- list(
    # gcFile = "results/tables/compiled_readGC_test.tsv",
    gcFile = "results/tables/compiled_readGC.tsv",
    # mdFile = "data/metadata.tsv",
    mdFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/results/starSummary/md_merged.tsv",
    outDir = "analysis/explore_GC"
)
io

# create outdir as needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}


### Format data
###############
# metadata
md <- read.table(io$mdFile, header = TRUE, row.names = NULL, sep = "\t")
md$rna_sampleID <- tstrsplit(md$sampleID, split = "_", fixed = TRUE)[[1]]

# gc distribution file
gc <- read.table(io$gcFile, header = TRUE, row.names = 1, sep = "\t")
colnames(gc) <- tstrsplit(colnames(gc), split = "_", fixed = TRUE)[[1]]
gc$pct_GC <- rownames(gc)
head(gc)
dim(gc)


########################################################
#####           Summarize GC Distributions         #####
########################################################

### Plot GC distributions
#########################
# melt table for ggplot
gc.melted <- reshape2::melt(gc, value.name = "count")
gc.melted$variable <- as.factor(gc.melted$variable)
gc.melted$pct_GC <- as.numeric(gc.melted$pct_GC)
gc.melted <- na.omit(gc.melted)
head(gc.melted)

# get the mean pct_GC for each sample
gcsMean <- list()
for (mysample in levels(gc.melted$variable)) {
    sub <- gc.melted[which(gc.melted$variable == mysample), ]
    sub$prod <- sub$pct_GC * sub$count
    gcsMean[[mysample]] <- ( sum(sub$prod) / sum(sub$count) )
}

# melt for plotting
gcsMean.melted <- reshape2::melt(as.data.frame(gcsMean))
colnames(gcsMean.melted) <- c("sampleID", "mean_pctGC")
m <- match(gcsMean.melted$sampleID, md$rna_sampleID)
gcsMean.melted$
head(gcsMean.melted)

# add mean percent GC to metadata
md$meanGC <- gcsMean.melted$value[match(gcsMean.melted$variable, md$rna_sampleID)]

# add other variables to melted table
m <- match(gc.melted$variable, md$rna_sampleID)
gc.melted$RINcat <- md$RINcat[m]
gc.melted$RINcat <- factor(gc.melted$RINcat, levels = c("low", "med",'high'))
gc.melted$RIN <- md$Novogene_RIN[m]
gc.melted$pctUnmapped <- md$pctUnmapped[m]
gc.melted$numMappedReads <- md$numUniqMappedReads[m]
head(gc.melted)

# histogram, colored by sample
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

# density plots
# colored by sample
ggplot(gc.melted, aes(x = pct_GC, color = variable)) +
    geom_density() +
    xlab("Percent GC") +
    ylab("Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 8), 
          axis.text.y = element_text(size = 8),
          legend.position = "None")
ggsave(paste(io$outDir, "gcDensity_bySample.png", sep = "/"), device = "png")

# density plots
# colored by RINcat
ggplot(gc.melted, aes(x = pct_GC, color = RINcat)) +
    geom_density() +
    xlab("Percent GC") +
    ylab("Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 8), 
          axis.text.y = element_text(size = 8),
          legend.position = "None")
ggsave(paste(io$outDir, "gcDensity_byRINcat.png", sep = "/"), device = "png")

# density plots
# colored by RIN
ggplot(gc.melted, aes(x = pct_GC, color = RIN)) +
    geom_density() +
    xlab("Percent GC") +
    ylab("Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 8), 
          axis.text.y = element_text(size = 8),
          legend.position = "None")
ggsave(paste(io$outDir, "gcDensity_byRINcat.png", sep = "/"), device = "png")

# density plots
# colored by %mapped
ggplot(gc.melted, aes(x = pct_GC, color = RINcat)) +
    geom_density() +
    xlab("Percent GC") +
    ylab("Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 8), 
          axis.text.y = element_text(size = 8),
          legend.position = "None")
ggsave(paste(io$outDir, "gcDensity_byRINcat.png", sep = "/"), device = "png")

# facet grid by sample
ggplot(gc.melted, aes(x = pct_GC, color = variable)) +
    geom_density() +
    facet_wrap(~variable) +
    geom_vline(aes(xintercept = value), gcsMean.melted, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "None")
ggsave(paste(io$outDir, "gcDensity_facet.png", sep = "/"), device = "png")


########################################################
#####             Correlate Mean %GC               #####
########################################################

### Pearson correlation test
############################
# apply Pearson correlation test
cor.test <- rcorr()