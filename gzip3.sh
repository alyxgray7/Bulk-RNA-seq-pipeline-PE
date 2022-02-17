#!/usr/bin/env bash
#SBATCH --partition=exacloud
#SBATCH --job-name=gzip3
#SBATCH --output=gzip3_%j.out
#SBATCH --error=gzip3_%j.err
#SBATCH --qos=long_jobs
#SBATCH --time=3-00:00:00

# variables
wDir=/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021

# gzip fasta sequence counts
ls -1 $wDir/samples/star/*/Unmapped.out.mate* | grep -v ".gz" | while read line ; do gzip $line ; done
