library(Glimma)
library(limma)
library(DESeq2)

project_id = snakemake@params[['project_id']]

rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))


out_path = file.path(getwd(),'results','diffexp')
if(!(file.exists( out_path ))) {
  print(paste("mkdir:", out_path))
  dir.create(out_path, FALSE, TRUE)  
}
# dir.create(out_path)

rds = readRDS(rds)
groups.df = as.data.frame(colData(rds))
glMDSPlot(rds, top = 1000, groups=groups.df,path=out_path,html=paste(project_id,'mds_plot',sep='.'),launch=FALSE)
