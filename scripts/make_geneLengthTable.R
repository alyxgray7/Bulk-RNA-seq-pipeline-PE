################################################################################
#####         Produce geneLength look up table from GTF reference file     #####
################################################################################

# Reads in a GTF file and calculates the gene lengths for each feature provided
# in the reference.
# Outputs a look up table with columns including ensemble_gene_id, external_gene_name,
# length_bp, and length_kb.
# Table is used in scripts/estSaturation.R

################################################################################
#####                            Set up parameters                         #####
################################################################################

### Command line arguments
##########################
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Calling the help function
help <- function() {
  cat("get_geneLengthsTable.R : Produces a gene lengths look up table, which is
                         used in estSaturation.R.
      
      - Outputs :
      
      1) gene lengths look up table at data/geneLengths.tsv.
      ")
  cat("\n")
  cat("Usage : \n")
  cat("--gtfFile          : GTF file from config['gtf_file'].                   [ required ]
      \n")
  cat("--annoFile         : Annotation file; from config['filter_anno'].        [ required ]
      \n")
  cat("--outDir           : Path to directory where results and data are saved. [ required ]
      \n")
  cat("\n")
  q()
}

# save values of each argument
if(!is.na(charmatch("--help", args)) || !is.na(charmatch("-h", args))){
  help()
} else {
  gtfFile <- sub('--gtfFile=', '', args[grep('--gtfFile=', args)])
  annoFile <- sub('--annoFile=', '', args[grep('--annoFile=', args)])
  outDir <- sub('--outDir=', '', args[grep('--outDir=', args)])
}
io <- list(
  gtfFile = gtfFile,
  annoFile = annoFile,
  outDir = outDir
)
io

# create outdir as needed
if(!(file.exists( outDir ))) {
  print(paste("mkdir:", outDir))
  dir.create(outDir, FALSE, TRUE)  
}


# ### Debugging on exa
# #############
# io <- list(
#   gtfFile = "/home/groups/CEDAR/anno/gtf/hg38_ens94.chr.gtf",
#   annoFile = "/home/groups/CEDAR/anno/biomaRt/hg38.Ens_94.biomaRt.geneAnno.Rdata",
#   outDir = "/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/pilot2/Bulk-RNA-seq-pipeline-PE/data"
# )
# io


### Libraries
#############
library(data.table)

################################################################################
#####                    Get gene lengths from GTF                         #####
################################################################################

### Run bash command to produce temporary table
###############################################
# command to run
a <- paste0("cat ", io$gtfFile, " ")
b <- " grep -v \'#' | grep \'exon' | awk -vOFS=\'\\t\' \'{L=$5-$4; a[$10]+=L}END{for(i in a){print i, a[i]}}' | tr -d \'\";'"
cmd <- paste(a, b, sep = "|")
cmd

# run command and save output as a datatable
tmp <- as.data.table(system(cmd, intern = TRUE))

# split table by the tab and make final changes
tmp[, c("ensembl_gene_id", "length_bp") := tstrsplit(V1, split = "\t", fixed = TRUE)]
final <- tmp[, c(2:3)]


### Add external gene names and save output
###########################################
# load annotation
load(io$annoFile, verbose = TRUE)

# translate ensembl_id to external_gene_name
m <- match(final$ensembl_gene_id, anno$ensembl_gene_id)
stopifnot(which(is.na(m)) == integer(0))
final$external_gene_name <- anno$external_gene_name[m]

# convert bp to kb
final$length_bp <- as.numeric(final$length_bp)
final$length_kb <- final$length_bp / 1000

# reorganize columns and save output
final <- final[, c(1,3,2,4)]
head(final)
write.table(final, paste(io$outDir, "geneLengths.tsv", sep = "/"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
