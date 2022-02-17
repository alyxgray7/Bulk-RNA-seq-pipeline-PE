#!/usr/bin/env bash
#SBATCH --partition=exacloud
#SBATCH --job-name=gzip1
#SBATCH --output=gzip1_%j.out
#SBATCH --error=gzip1_%j.err
#SBATCH --qos=long_jobs
#SBATCH --time=1-00:00:00

# variables
wDir=/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021

# gzip overrepresented sequence counts
ls -1 $wDir/data/unmappedSeqs/*.txt | while read line ; do gzip $line ; done
