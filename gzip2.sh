#!/usr/bin/env bash
#SBATCH --partition=exacloud
#SBATCH --job-name=gzip2
#SBATCH --output=gzip2_%j.out
#SBATCH --error=gzip2_%j.err
#SBATCH --qos=long_jobs
#SBATCH --time=1-00:00:00

# variables
wDir=/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021

# gzip fasta sequence counts
ls -1 $wDir/samples/star/*/*.fa | while read line ; do gzip $line ; done
