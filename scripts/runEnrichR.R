#############################################################
#####       Script to run EnrichR pathway enrichment    #####
#############################################################

### Set up
##########
# command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
help <- function(){
  cat("RNAseq_readQC.R :
  
  - Outputs :
  
  1) Gene enrichment tables and plots for included databases. \n
  2) Databases used: \n
        - KEGG_2021_Human \n
        - GO_Biological_Process_2021 \n
        - GO_Cellular_Component_2021 \n
        - GO_Molecular_Function_2021 \n
      ")
  cat("\n")
  cat("Usage : \n")
  cat("--degFile        : [ required ] Path to DEG table output by pairwise comparisons.
    \n")
  cat("--metaFile       : [ required ] Path to config['omic_meta_data']
    \n")
  cat("--annoFile       : [ required ] Path to config['filter_anno']
    \n")
  cat("--outDir         : [ required ] Path to the directory to save the figures in.
    \n")
  cat("--sampleID       : [ required ] Column specifying the sample IDs. Specified by config['sample_id'].
      \n")
  cat("--padj           : [ required ] Desired FDR cutoff.
  \n")
  cat("--FC             : [ required ] Desired fold change cutoff (not log2!)\n")
  cat("\n")
  q()
}

# save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  degFile <- sub('--degFile=', '', args[grep('--degFile=', args)])
  metaFile <- sub('--metaFile=', '', args[grep('--metaFile=', args)])
  annoFile <- sub('--annoFile=', '', args[grep('--annoFile=', args)])
  outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
  sampleID <- sub('--sampleID=', '', args[grep('--sampleID=', args)])
  padj <- sub('--padj=', '', args[grep('--padj=', args)])
  FC <- sub('--FC=', '', args[grep('--FC=', args)])
}

# std in/out
io <- list(
    degFile = degFile
    , metaFile = metaFile
    , annoFile = annoFile
    , outDir = outDir
    , sampleID = sampleID
    , padj = as.numeric(padj)
    , FC = as.numeric(FC)
)

# for debugging on exa
# io <- list(
#     degFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_02162022/results/diffexp/pairwise/1_Case-vs-4_ScreenNegative.diffexp.tsv"
#     , metaFile = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_02162022/data/pltRNAseq_metadata_02162022.tsv"
#     , annoFile = "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_94.biomaRt.geneAnno.Rdata"
#     , outDir = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_02162022/results/diffexp/pairwise/enrichR_test"
#     , sampleID = "rnaSampleID"
#     , padj = 0.05
#     , FC = 1.5
# )

# get contrast from filename
contrast <- tail(strsplit(io$degFile, split = "/", fixed = TRUE)[[1]], 1)
io$contrast <- gsub(".diffexp.tsv", "", contrast)
io

# create outdir as needed
if(!(file.exists( io$outDir ))) {
  print(paste("mkdir:", io$outDir))
  dir.create(io$outDir, FALSE, TRUE)  
}

# libraries
library(ggplot2)
library(enrichR)
library(data.table)


### Format data
###############
# anno
load(io$annoFile, verbose = TRUE)

# metadata
md <- read.table(io$metaFile, header = TRUE, row.names = NULL, sep = "\t")
rownames(md) <- md[[io$sampleID]]
head(md)

# table including DEGs from pairwise contrast
deg <- read.table(io$degFile, header = TRUE, row.names = 1, sep = "\t")

# determine how genes are labeled
annoType <- c()
if (length(grep("ENSG", rownames(deg))) == 0) {
  print("counts are external gene ids")
  annoType <- "external_gene_name"
  deg[["external_gene_name"]] <- rownames(deg)
  rownames(deg) <- NULL
} else {
  print("counts are ensembl ids")
  annoType <- "ensembl_gene_id"
  m <- match(rownames(deg), anno[[annoType]])
  deg[["external_gene_name"]] <- anno[["external_gene_name"]][m]
  deg[["ensembl_gene_id"]] <- rownames(deg)
  rownames(deg) <- NULL
}
head(deg)

# get list of up and down DEGs
up <- which(deg$padj < io$padj & deg$log2FoldChange >= log2(io$FC))
up <- deg[up,]
down <- which(deg$padj < io$padj & deg$log2FoldChange <= -log2(io$FC))
down <- deg[down,]
geneList <- list(
    up = unique(up[["external_gene_name"]])
    , down = unique(down[["external_gene_name"]])
)


### EnrichR
###########
setEnrichrSite("Enrichr")

# grab databases
dbs2get <- c("KEGG_2021_Human"
             , "GO_Biological_Process_2021"
             , "GO_Cellular_Component_2021"
             , "GO_Molecular_Function_2021")

# loop through up and down genes
for (myList in names(geneList)) {
    print(paste0("Running EnrichR for ", myList, " genes"))
    enriched <- enrichr(geneList[[myList]], dbs2get)

    # loop through each database search
    for (db in names(enriched)) {

        # save results
        write.table(enriched[[db]], paste0(io$outDir, "/", io$contrast, "-", db, ".", myList, "FC.", io$FC, ".adjp.", io$padj, ".tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
        plotEnrich(enriched[[db]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
        ggsave(paste0(io$outDir, "/", io$contrast, "-", db, ".", myList, "FC.", io$FC, ".adjp.", io$padj, ".pdf"), device = "pdf")
    }
}
