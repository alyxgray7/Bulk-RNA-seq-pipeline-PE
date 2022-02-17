#!/usr/bin/env bash
#SBATCH --partition=exacloud
#SBATCH --job-name=test_geneBody
#SBATCH --output=test_geneBody_%j.log
#SBATCH --error=test_geneBody_%j.log
#SBATCH --qos=long_jobs
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=10Gb
#SBATCH --cpus-per-task=2


# load environment
source /home/groups/CEDAR/grayaly/miniconda3/etc/profile.d/conda.sh
conda activate rseqc2
echo "activated environment"

# set up vars
wDir=/home/groups/CEDAR/grayaly/projects/platelet/plt-rnaseq/full-cohort/Bulk-RNA-seq-pipeline-PE_12092021

# run command
geneBody_coverage.py -r /home/groups/CEDAR/anno/gtf/hg38_ens94.chr.bed \
-i $wDir/samples/star/F4_InSitu_bam/Aligned.sortedByCoord.out.bam \
-o $wDir/test_geneBody
