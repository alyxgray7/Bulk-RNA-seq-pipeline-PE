library(Glimma)
library(limma)
library(DESeq2)
library(stringr)

condition = snakemake@params[['condition']]
cat(sprintf(c('Condition: ',condition,'\n')))

ma_plot_ = snakemake@output[['ma_plot']]
volcano_plot_ = snakemake@output[['volcano_plot']]
ma_plot = str_sub(tail(strsplit(ma_plot_,'/')[[1]],n=1),1,-6)
volcano_plot = str_sub(tail(strsplit(volcano_plot_,'/')[[1]],n=1),1,-6)


contrast = c(condition, snakemake@params[["contrast"]])
rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

# create outdir as needed
out_path = file.path(getwd(),'results','diffexp', 'covariate')
if(!(file.exists( out_path ))) {
  print(paste("mkdir:", out_path))
  dir.create(out_path, FALSE, TRUE)  
}
# dir.create(out_path)
print(out_path)

rds = readRDS(rds)
groups.df = as.data.frame(colData(rds))


#### by contrasts
res <- results(rds, contrast=contrast)
res$padj[is.na(res$padj)] = 1

rnaseq = as.data.frame(counts(rds, normalized=T))
genes = as.data.frame(row.names(res))
colnames(genes) = 'GeneID'

status_frame = res[,c('log2FoldChange','padj')]
status_frame['status'] = 0
status_frame$padj[is.na(status_frame$padj)] = 1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange < 0 ,'status'] = -1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange > 0 ,'status'] = 1

glMDPlot(res, anno=genes, status=status_frame$status, samples=colnames(rnaseq), 
         counts=log2(rnaseq + 0.0001),
        groups=groups.df[[condition]], main=strsplit(res@elementMetadata$description[2],': ')[[1]][2],
         transform=F, side.ylab='Log2-expression',launch=FALSE,side.main='GeneID', html = ma_plot, path=out_path)

## Volcano plot
glXYPlot(x=res$log2FoldChange, y=-log10(res$pvalue), xlab="logFC", ylab="logodds",path=out_path,
         status=status_frame$status, launch=FALSE,counts=log2(rnaseq + 0.0001), groups=groups.df[[condition]],
         anno=genes, main=strsplit(res@elementMetadata$description[2],': ')[[1]][2], html = volcano_plot)

