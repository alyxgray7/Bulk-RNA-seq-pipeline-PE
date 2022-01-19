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
library(EDASeq)
# library(Hmisc)
# library(pheatmap)
# library(corrplot)


### Std in/out
##############
io <- list(
    # gcFile = "results/tables/compiled_readGC_test.tsv",
    # gcFile = "results/tables/compiled_readGC.tsv",
    # mdFile = "data/metadata.tsv",
    mdFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/results/starSummary/md_merged.tsv",
    annoFile = "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_94.biomaRt.geneAnno.Rdata",
    gcFile = "/home/groups/CEDAR/grayaly/annotations/genomes/hg38/hg38.gc_content.tsv",
    countsFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/data/platelet_full-cohort_counts.txt",
    outDir = "/home/groups/CEDAR/grayaly/projects/platelet/analysis/edaseq",
    sampleID = "sampleID"
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
head(md[[io$sampleID]])

# raw counts
# raw <- read.table(io$countsFile, header = TRUE, row.names = NULL, sep = "\t")
# genesAll <- raw$Genes
# raw <- raw[,-1]
# stopifnot(colnames(raw) %in% md[[io$sampleID]])
# dim(raw)
raw <- read.table(io$countsFile, header = TRUE, row.names = 1, sep = "\t")
raw <- raw[, which(colnames(raw) %in% md[[io$sampleID]])]
raw <- raw[which(apply(raw, 1, mean) <= 5), ]
dim(raw)

# # only consider genes with mean counts above 5
# genes2keep <- genesAll[which(apply(raw, 1, mean) >=5)]
# raw <- raw[which(apply(raw, 1, mean) >= 5), ]
# stopifnot(length(genes2keep) == nrow(raw))
# dim(raw)

# annotations
anno <- get(load(io$annoFile))
if (length(grep("ENSG", rownames(raw))) == 0) {
  print("counts are external gene ids")
  annoType <- "external_gene_name"
} else {
  print("counts are ensembl ids")
  annoType <- "ensembl_gene_id"
}
annoType

# gc percentages
gc <- read.table(io$gcFile, header = TRUE, row.names = NULL, sep = "\t")
gc$gene_length <- gc$gene_end_bp - gc$gene_start_bp
cols2keep <- c(annoType, "percent_GC", "gene_length")
gc <- gc[, which(colnames(gc) %in% cols2keep)]
gc <- gc[which(!gc[[annoType]] == ""), ]

# subset for genes included in raw counts
gc.sub <- gc[which(gc[[annoType]] %in% rownames(raw))]
#gc.sub <- gc[which(gc[[annoType]] %in% genes2keep), ]
dim(gc.sub)


### Start exploratory data analysis (EDAseq)
############################################
# ensure that the sizes between feature data and raw counts match


# create SeqExpressionSet object
data <- newSeqExpressionSet(
  counts = as.matrix(raw),
  featureData = gc.sub,
  phenoData = md
)
feature <- data.frame(gc =)


