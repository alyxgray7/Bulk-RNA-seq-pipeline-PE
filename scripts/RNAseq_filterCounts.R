args <- commandArgs()
print(args)
help <- function() {
    cat("RNAseq_filterCounts.R :
  
  - Outputs :
  
  1) Tab separated filtered counts table for biotypes specified in config['biotypes'] and removal of mitochondrial genes as specified by config['mito']. 
            Output: {project_id}_{counts_file_name}.filt.txt
  2) If included, a second counts table with ERCC genes is printed to as a tab separated table. Process is specified by config['ERCC'].
            Output: {project_id}_{counts_file_name}.ercc.txt
      ")
  cat("\n")
  cat("Usage : \n")
  cat("--annoFile     : [ required ] Path to config['filter_anno']
      \n")
  cat("--countsFile   : [ required ] Path to unfiltered counts matrix.
      \n")
  cat("--biotypes     : [ required ] Comma separated biotypes to include in filtered counts table.")
  cat("--mito         : [ required ] Binary way to signal if mitochondrial genes should be included in filtered counts table (0) or removed (1).")
  cat("--ercc         : [ required ] Binary way to signal if ERCC genes were included (1) in counts table or not (0).")
  cat("--outDir       : [ required ] Desired directory to write table outputs.")
  cat("\n")
  q()
}

# save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  annoFile <- sub('--annoFile=', '', args[grep('--annoFile=', args)])
  countsFile <- sub('--countsFile=', '', args[grep('--countsFile=', args)])
  biotypes <- sub('--biotypes=', '', args[grep('--biotypes=', args)])
  mito <- sub('--mito=', '', args[grep('--mito=', args)])
  ercc <- sub('--ercc=', '', args[grep('--ercc=', args)])
  outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
}

# check input files in logs/.out file
io <- list(
  annoFile = annoFile,
  countsFile = countsFile,
  biotypes = biotypes,
  mito = mito,
  ercc = ercc,
  outDir = outDir
)
io

# create outdir as needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}

# annoFile = snakemake@params[['anno']]
# biotypes <- snakemake@params[['biotypes']]
# countsFile <- snakemake@input[['countsFile']]
# mito <- snakemake@params[['mito']]

# debug on exa
# io <- list(
#     annoFile = "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_94.biomaRt.geneAnno.Rdata"
#     , biotypes  = ""
#     , countsFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/data/platelet_full-cohort_genecounts.txt"
#     # , countsFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/data/platelet_full-cohort_counts.txt"
#     , mito = "1"
#     , ercc = "1"
#     , outDir = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021/data"
# )
# io


### Load data
#############
# counts
print("Loading counts table")
print(io$countsFile)

# check if an rda file or tab sep
if(grepl('rda|RData|Rdata', io$countsFile)){
    counts <- get(load(file = io$countsFile))
    head(counts[,1:5])
    dim(counts)
}
if(grepl('txt|tsv', io$countsFile)){
    counts <- read.delim(file = io$countsFile)
    head(counts[,1:5])
    dim(counts)
}

# anno file
print("Loading annotation table")
print(io$annoFile)

# check if an rda file or tab sep
if(grepl('rda|RData|Rdata',annoFile)){
    anno <- get(load(file = io$annoFile))
}
if(grepl('txt|tsv',annoFile)){
    anno <- read.delim(file = io$annoFile)
}

# check gene annotation type in counts table
annoType <- c()
if (length(grep("ENSG", counts[,1])) == 0) {
  print("counts are external gene ids")
  annoType <- "external_gene_name"
} else {
  print("counts are ensembl ids")
  annoType <- "ensembl_gene_id"
}
annoType

### Filter
##########
# filter annotations for biotypes we're interested in keeping
if (length(strsplit(io$biotypes, split = "\\,")[[1]]) > 0) {
    anno.sub <- anno[which(anno$gene_biotype %in% strsplit(io$biotypes, split = "\\,")[[1]]), ]
    counts.sub <- counts[which(counts[,1] %in% unique(paste(anno.sub[[annoType]]))), ]
} else {
    print("no biotypes provided")
    counts.sub <- counts
}

# filter mitochondrial genes if desired
if (io$mito == 1 & annoType == "external_gene_name") {
    print("tossing MT- genes")
    counts.sub <- counts.sub[grep("^MT-", paste(counts.sub$Genes), invert=TRUE), ]
}
if (mito == 1 & annoType != "external_gene_name") {
    print("gene names are ensembl_ids; tossing MT- genes")
    mtgenes <- anno[[annoType]][grep("MT-", anno$external_gene_name)]
    counts.sub <- counts.sub[which(!counts.sub[,1] %in% mtgenes), ]
}

# save filtered counts table
filename <- tail(strsplit(io$countsFile, split = "/", fixed = TRUE)[[1]], 1)
filename <- sub(".txt", ".filt.txt", filename)
write.table(counts.sub, file=paste0(io$outDir, "/", filename), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# filter ERCC genes if desired
if (io$ercc == 1) {
    print("filtering ERCC genes")
    counts.ercc <- counts[grep("ERCC", counts[,1]), ]

    # save output
    filename <- tail(strsplit(io$countsFile, split = "/", fixed = TRUE)[[1]], 1)
    filename <- sub(".txt", ".ercc.txt", filename)
    write.table(counts.ercc, file = paste0(io$outDir, "/", filename), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

print("Finished script")

