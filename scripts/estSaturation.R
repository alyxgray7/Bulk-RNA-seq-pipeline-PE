################################################################################
#####         Estimate sequencing saturation from raw counts matrix        #####
################################################################################

# Estimates the sequencing saturation and for a given RNA-seq experiment.
# Produces library complexity curve plots.
# Calls functions `estimate_saturation()` and `plot_saturation_curves()` 
# from RNAseQC.
# (https://rdrr.io/github/mjdufort/RNAseQC/man/estimate_saturation.html)


################################################################################
#####                            Set up parameters                         #####
################################################################################

### Command line arguments
##########################
# set up args
args <- commandArgs(trailingOnly = TRUE)
print(args)

# calling the help function
help <- function() {
  cat("estSaturation.R : Estimates the sequencing saturation for a given RNA-seq experiment.
                         Uses estimate_saturation and plot_saturation_curves from RNAseQC.
                         (https://rdrr.io/github/mjdufort/RNAseQC/man/estimate_saturation.html)
                         Plots the number of features identified for each sample as well as the
                         distribution of the number of features by contrast group.
      
      - Outputs :
      
      1) Plot of sequencing saturation curves, colored by sample.
      2) Plot of sequencing saturation curves, split by contrast group.
      3) Plot of sequencing saturation curves, split by sample.
      4) Plot of the number of features per sample with standard error bars, colored by contrast group.
      5) Plot of the distribution of feature number by contrast group.
      6) Plot of the distribution of sequencing saturation variance by contrast group.
      7) Tab-separated table of values returned from RNAseQC::estimate_saturation().
      8) Tab-separated table of sequencing summary statistics.
      9) Tab-separated table of raw counts normalized to FPKM.
      ")
  cat("\n")
  cat("Usage : \n")
  cat("--countsFile       : Counts file produced by the pipeline.               [ required ]
      \n")
  cat("--mdFile           : Metadata file; from config['omic_meta_data'].       [ required ]
      \n")
  cat("--geneLengthsFile  : table of gene lengths produced from pipeline        [ required ]
                            step make_geneLengthsTable.
      \n")
  cat("--outDir           : Path to directory where results and data are saved. [ required ]
      \n")
  cat("--minCount         : Sets the threshold for the minimum count needed for [ required ]
                            a feature to be 'expressed'. Default = 1.
      \n")
  cat("--contrast     : The column name in your config['omic_meta_data'] file,  [ required ]
                        this is the characteristic you would like to do DE on.
                        Example: diagnosis, geotype, etc. (used to color the plots by.)
      \n")
  cat("\n")
  q()
}

# save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  countsFile <- sub('--countsFile=', '', args[grep('--countsFile=', args)])
  mdFile <- sub('--mdFile=', '', args[grep('--mdFile=', args)])
  geneLengthsFile <- sub('--geneLengthsFile=', '', args[grep('--geneLengthsFile=', args)])
  outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
  minCount <- as.numeric(sub('--minCount=', '', args[grep('--minCount=', args)]))
  contrast <- sub('--contrast=', '', args[grep('--contrast=', args)])
}

# create outdir as needed
if(!(file.exists( outDir ))) {
  print(paste("mkdir:", outDir))
  dir.create(outDir, FALSE, TRUE)  
}

# save std in/out
io <- list(
  countsFile = countsFile,
  mdFile = mdFile,
  geneLengthsFile = geneLengthsFile,
  outDir = outDir,
  minCount = as.numeric(minCount),
  contrast = contrast
)
io

### Debugging
# local
# io <- list(
#   countsFile = "~/Box/PLATELET RNA & MULTIOMICS/DATA/Pilots_combo/combo_counts.txt",
#   mdFile = "~/Box/PLATELET RNA & MULTIOMICS/DATA/Pilots_combo/combo_metadata.txt",
#   geneLengthsFile = "~/CEDAR/annotations/hg38_ens94.chr_transcriptLengths.txt",
#   outDir = "~/CEDAR/projects/platelet/pancreas/rnaseq/analysis/pilots_combo/estSaturation",
#   minCount = 2,
#   contrast = "RINcat"
# )

# # exa
# io <- list(
#   countsFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/pilot2/Bulk-RNA-seq-pipeline-PE/data/platelet_pilot2_counts.txt",
#   mdFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/pilot2/Bulk-RNA-seq-pipeline-PE/data/metadata.tsv",
#   geneLengthsFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/pilot2/Bulk-RNA-seq-pipeline-PE/data/geneLengths.tsv",
#   outDir = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/pilot2/Bulk-RNA-seq-pipeline-PE/results/estSaturation",
#   minCount = 2,
#   contrast = "RINcat"
# )
# io


### Libraries
#############
# all necessary packages
packages <- c("ggplot2", "ggthemes", "scales", "RColorBrewer", "remotes",
              "biomaRt", "edgeR", "ComplexHeatmap", "countSubsetNorm", "RNAseQC")

# load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) { 
    if (!require(x, character.only = TRUE)) {
      if (x == "RNAseQC" || x == "countSubsetNorm") {
        # requires remote install with github
        remotes::install_github(paste("mjdufort", x, sep = "/"), 
                                upgrade = c("never"))
        library(x, character.only = TRUE)
      } else {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  }
)

################################################################################
#####                                Functions                             #####
################################################################################

### Calculate fpkm from raw counts
##################################
# function to calculate rpkm/fpkm
# matrix = counts matrix
# df = geneLengths look up table
fpkm.fx <- function(matrix, df, annoType) {
  
  # normalize per million fragments
  pm <- apply(matrix, 2, sum) / 10^6
  
  # divide counts by per million scaling factor
  cpm.mat <- matrix / pm
  
  # ensure the correct order based on annotation type
  m <- match(rownames(matrix), df[, annoType])
  
  # divide counts per million by length of gene in kilobases
  rpkm.mat <- cpm.mat / df$length_kb[m]
  
  return(rpkm.mat)
}

################################################################################
#####                              Format data                             #####
################################################################################

### Raw counts
##############
# read in raw counts matrix
counts <- read.table(io$countsFile, header = TRUE, row.names = NULL, sep = "\t")

# reformat into matrix and remove uninformative features
counts.mat <- as.matrix(counts[, -1])
rownames(counts.mat) <- counts[, 1]
counts.mat <- counts.mat[which(rowSums(counts.mat) >= 1), ]

# check how gene names are annotated in counts
annoType <- c()
if (length(grep("ENSG", rownames(counts.mat))) == 0) {
  print("counts are external gene ids")
  annoType <- "external_gene_name"
} else {
  print("counts are ensembl ids")
  annoType <- "ensembl_gene_id"
  rownames(counts.mat) <- sub("\\..*$", "", rownames(raw)) # be sure ens Ids are unique
}


### GeneLengths lookup table
############################
# read in geneLengths file and rename columnns
geneLengths <- read.table(io$geneLengthsFile, header = TRUE, row.names = NULL, sep = "\t")
head(geneLengths)


### Normalize raw counts per kilobase million
#############################################
# PE data uses "fragments" instead of "reads"
fpkm.mat <- fpkm.fx(counts.mat, geneLengths, annoType)
print("Normalized (FPKM) values")
head(fpkm.mat)

# sanity check -- needs more ambiguous solution
if (annoType == "external_gene_name") {
  c <- counts.mat["DDX11L1", 1]
  pm <- (sum(counts.mat[,1])/10^6)
  cpm <- c/pm
  kb <- geneLengths[grep("DDX11L1", geneLengths[, annoType]), ]$length_kb[2]
  stopifnot(fpkm.mat["DDX11L1", 1] == cpm / kb)
}

# save for future use
write.table(fpkm.mat, paste(io$outDir, "counts_fpkm.tsv", sep = "/"), 
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


### Meta data
#############
# read in metadata file
md <- read.table(io$mdFile, header = TRUE, row.names = NULL, sep = "\t")
head(md)


### Get sequencing statistics
#############################
# initialize dataframe
stats <- data.frame(
  SampleID = md[,1]
)

# add counts summaries
stats$sum_fpkm <- round(apply(fpkm.mat, 2, sum), 3)
stats$mean_fpkm <- round(apply(fpkm.mat, 2, mean), 3)
stats$med_fpkm <- round(apply(fpkm.mat, 2, median), 3)
stats$sd_fpkm <- round(apply(fpkm.mat, 2, sd), 3)
stats$cv_fpkm <- round((apply(fpkm.mat, 2, sd) / apply(fpkm.mat, 2, mean))*100, 3)
stats$var_fpkm <- round(apply(fpkm.mat, 2, var), 3)
stats

# save table
write.table(stats, paste(io$outDir, "summaryStatistics.tsv", sep = "/"), 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


################################################################################
#####                           Estimate saturation                        #####
################################################################################

### Estimate sequencing saturation 
#####################################
# call RNASeQC function on normalized counts
# sample nFeatures at 20 different points between 0 and the maximum nFeatures from all samples
# repeat 100 times
print("Begin estimating saturation")
satEst <- estimate_saturation(
  counts = fpkm.mat,
  method = "sampling",
  ndepths = 20,
  nreps = 100,
  min_counts = io$minCount,
  verbose = TRUE
)


### Reformat results for saving
###############################
# add contrast information
m <- match(satEst$sample, md[, 1])
satEst$contrast <- md[, io$contrast][m]

# save results
tosave <- satEst
colnames(tosave) <- c("SampleID", "readDepth", "nFeatures", "variance", io$contrast)
write.table(tosave, paste(io$outDir, "estSaturation_results.tsv", sep = "/"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
print("Finished estimating saturation")


################################################################################
#####                             Summary plots                            #####
################################################################################

### Plotting parameters
#######################
# number of colors to create in palette
colCount <- length(unique(satEst$sample))

# set palette
myPalette <- "Set2"


### Plot saturation curves
##########################
# include all samples and color by sample
plot_saturation_curve(satEst, plot_terminal_points = FALSE) +
  ggtitle("Sequencing Saturation Curve") +
  ylab("Number of features detected") +
  xlab("Number of mapped reads") +
  scale_x_continuous(labels = scientific) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, myPalette))(colCount)) +
  theme_bw()
ggsave(paste(io$outDir, "saturationCurve_bySample.png", sep = "/"), device = "png")

# split by sample
plot_saturation_curve(satEst) +
  ggtitle("Sequencing Saturation Curve") +
  ylab("Number of features detected") +
  xlab("Number of mapped reads") +
  scale_x_continuous(labels = scientific) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, myPalette))(colCount)) +
  theme_bw() +
  facet_wrap( ~ sample) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste(io$outDir, "facet_saturationCurve_bySample.png", sep = "/"), device = "png")

# split by contrast
plot_saturation_curve(satEst, 
                      design = satEst, 
                      design_id_col = c("contrast"), 
                      color_lines_by_var = c("contrast"),
                      plot_terminal_points = FALSE) +
  facet_grid(~ contrast) +
  ggtitle("Sequencing Saturation Curve") +
  ylab("Number of features detected") +
  xlab("Number of mapped reads") +
  scale_x_continuous(labels = scientific) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, myPalette))(colCount)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste(io$outDir, "facet_saturationCurve_byContrast.png", sep = "/"), device = "png")


### Plot variance distributions
###############################
# violin plot, colored by contrast
ggplot(satEst, aes(x = contrast, y = sat.var, fill = contrast)) +
  #facet_grid(~ contrast)
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white") +
  ggtitle("Distribution of sequencing saturation variance") +
  xlab("Contrast") +
  ylab("Variance") +
  theme(legend.position = "none") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, myPalette))(colCount)) +
  theme_bw()
ggsave(paste(io$outDir, "violin_saturationVariance_byContrast.png", sep = "/"), device = "png")


### Plot number of features
###########################
# create dataframe for plotting
toplot <- data.frame(
  sample = unique(satEst$sample)
)

# add number of features and its variance for each sample
nFeatures <- c()
variance <- c()
for (i in toplot$sample) {
  
  # get maximum features identified at the max depth
  value <- max(satEst[which(satEst$sample == i), ]$sat, na.rm = TRUE)
  nFeatures <- append(nFeatures, value, length(nFeatures))
  
  # save position of max features
  index <- which(satEst$sat == value)
  
  # get calculated variance from the estimated max feature
  variance <- append(variance, satEst[index, ]$sat.var, length(variance))
}

# add vectors to dataframe
toplot$nFeatures <- nFeatures
toplot$variance <- variance

# add contrast information
m <- match(toplot$sample, md[, 1])
toplot$contrast <- md[, io$contrast][m]

# get the SD from variance
toplot$sd <- sqrt(toplot$variance)

# order table by number of features and contrast
toplot <- toplot[order(toplot$contrast, toplot$nFeatures, decreasing = TRUE), ]
toplot$sample <- factor(toplot$sample, levels = toplot$sample)
toplot

# barplot with SE bars
ggplot(toplot, aes(x = sample, y = nFeatures, fill = contrast)) +
  geom_col() +
  geom_errorbar(aes(ymin = nFeatures-sd, ymax = nFeatures+sd), 
                width = 0.2, position = position_dodge(0.9)) +
  ggtitle(paste0("Number of features (min_count = ", io$minCount, ")")) +
  xlab("Sample") +
  ylab("Number of features") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, myPalette))(colCount)) +
  theme_bw()
ggsave(paste(io$outDir, "barplot_nFeatures_bySample.png", sep = "/"), device = "png")

# violin plot of feature numbers by contrast group
ggplot(toplot, aes(x = contrast, y = nFeatures, fill = contrast)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white") +
  ggtitle("Distribution of feature number by contrast") +
  xlab("Contrast") +
  ylab("Number of features") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, myPalette))(colCount)) +
  theme(legend.position = "none") +
  theme_bw()
ggsave(paste(io$outDir, "violin_nFeatureDistribution_byContrast.png", sep = "/"), device = "png")

print("Finished script")
