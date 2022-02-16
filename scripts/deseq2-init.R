### Set up
##########
# command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

# help function
help <- function() {
  cat("deseq2-init.R : Initializes the deseq2 object.
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
  cat("--contrast         : [ required ] Sample group contrasts to build DE model.
  \n")
  cat("--threads          : [ required ] Number of threads for parallel processing. (Defaults to 6 unless the number of samples is > 100.)
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
  contrast <- sub('--contrast=', '', args[grep('--contrast=', args)])
  threads <- sub('--threads=', '', args[grep('--threads=', args)])
}

# create the contrast vector
contrast <- as.vector(strsplit(contrast, split = " ", fixed = TRUE)[[1]])

# save std in/out
io <- list(
  countsFile = countsFile,
  outDir = outDir,
  metaFile = metaFile,
  sampleID = sampleID,
  Type = linear_model,
  contrast = contrast,
  threads = threads
)
io

# debug on exa
# io <- list(
#   countsFile = "data/platelet_full-cohort_genecounts.filt.txt",
#   outDir = "results/diffexp",
#   metaFile = "data/md_merged_noI9.tsv",
#   sampleID = "rnaSampleID",
#   Type = "Group",
#   contrast = contrast,
#   threads = 
# )
# io

# set up vars
target <- io$contrast[1]
print(paste0("Target: ", target))

baseline <- io$contrast[2]
print(paste0("Baseline: ", baseline))

output <- paste0(io$outDir, "/", paste0(target, "-vs-", baseline), "_all.rds")
print(output)

rld_out <- paste0(io$outDir, "/", paste0(target, "-vs-", baseline), "_rlog_dds.rds")
print(rld_out)

# create outdir as needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}

# load libraries
library("dplyr")
library("DESeq2")
parallel <- FALSE
# if parallel
if (io$threads > 1) {
    library("BiocParallel")
    register(MulticoreParam(io$threads))
    parallel <- TRUE
}


### Load & format data
######################
# read in metadata table and order according to sampleID
md <- read.delim(file = io$metaFile, sep = "\t", stringsAsFactors = FALSE)
md <- md[order(md[[io$sampleID]]),]
head(md)

# read in counts table
subdata <- read.table(io$countsFile, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
subdata <- subdata[ , order(colnames(subdata))]
dim(subdata)
head(subdata[,1:5])

# extract only the Types that we want in further analysis
# keep only the sampleIDs and status informative columns
md <- md[which(
    md[[io$Type]] == baseline |
    md[[io$Type]] == target &
    md[[io$sampleID]] %in% colnames(subdata)
),]
head(md)
dim(md)

# keep only the sampleIDs of the types we have chosen in the metadata table above
rownames(md) <- md[[io$sampleID]]
md[[io$sampleID]] <- NULL
keep <- colnames(subdata)[colnames(subdata) %in% rownames(md)]
subdata <- subdata[, keep]
dim(subdata)

# check
stopifnot(rownames(md) == colnames(subdata))


### DESeq2
##########
# obtain the number of genes that meet padj<0.01 for reference line in histogram
dds <- DESeqDataSetFromMatrix(countData = subdata,
                              colData = md,
                              design = as.formula(paste('~' , io$Type)))

dds <- estimateSizeFactors(dds)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]

# normalization and pre-processing
dds <- DESeq(dds, parallel=parallel)
saveRDS(dds, file=output)

# obtain normalized counts
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file=rld_out)

print("Finished script.")