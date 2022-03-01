__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""


import datetime
import sys
import os
import pandas as pd
import json


timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

# configfile:"omic_config.yaml"
configfile:"omic_config_byRisk.yaml"
project_id = config["project_id"]

### RunALL
SAMPLES, = glob_wildcards("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam")
# SAMPLES, = glob_wildcards("samples/htseq_count/{sample}_genecount.txt")

### Testing
# SAMPLES = ['A10', 'A1', "A2", "A3"]

ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['R1','R2']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']
readQC_plotNames = ['biotype_barplot', 'geneAttributes_barplot', 'totalReads_arranged']
estSat_plotNames = ['barplot_nFeatures_bySample','facet_saturationCurve_bySample', 'facet_saturationCurve_byContrast','violin_nFeatureDistribution_byContrast', 'violin_saturationVariance_byContrast']
estSat_tableNames = ['counts_fpkm', 'estSaturation_results', 'summaryStatistics']


with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
rule_dirs.pop(rule_dirs.index('__default__'))


for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


# def format_plot_columns():
#     factors = config['meta_columns_to_plot'].keys()
#     reformat_factors = '"' + '","'.join(factors) + '"'
#     return 'c({})'.format(reformat_factors)

def format_plot_columns(wildcards=None):
    factors = config['meta_columns_to_plot'].keys()
    reformat_factors = ",".join(factors)
    return '"{}"'.format(reformat_factors)


def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(config["omic_meta_data"]) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")


rule all:
    input:
        # expand("samples/fastqc/{sample}/{sample}_{fastq_ext}_t.good_fastqc.zip", sample = SAMPLES, fastq_ext = fastq_ext),
        # expand("samples/fastqscreen/{sample}/{sample}_{fastq_ext}_t.good_screen.{fastqscreen_ext}", sample=SAMPLES, fastq_ext=fastq_ext, fastqscreen_ext=fastqscreen_ext),
        # expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        # expand("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        # expand("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai", sample = SAMPLES),
        # expand("samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam.bai", sample = SAMPLES),
        # expand("samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam", sample = SAMPLES),
        expand("samples/htseq_count/{sample}_genecount.txt", sample = SAMPLES),
        "data/{project_id}_genecounts.txt".format(project_id=config["project_id"]),
        # "data/{project_id}_genecounts_w_stats.txt".format(project_id=config["project_id"]),
        # "data/{project_id}_genecounts.filt.txt".format(project_id=config["project_id"]),
        # "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"]),
        # expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        # expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
        # expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        # expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        # "results/tables/read_coverage.txt",
        # expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        # expand("rseqc/geneBody_coverage/{sample}/{sample}.geneBodyCoverage.curves.pdf", sample = SAMPLES),
        # expand("results/readQC/{plot}.png", plot = readQC_plotNames),
        expand("results/diffexp/pairwise/{contrast}_all.rds", contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/pairwise/{contrast}_rlog_dds.rds", contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/pairwise/{contrast}.pca_plot.pdf", contrast = config["diffexp"]["contrasts"]),
        "results/diffexp/group/LRT_pca.pdf",
        "results/diffexp/group/MDS_table.txt",
        "results/diffexp/group/LRT_density_plot.pdf",
        expand(["results/diffexp/pairwise/{contrast}.qplot.pdf","results/diffexp/pairwise/{contrast}.qhist.pdf","results/diffexp/pairwise/{contrast}.qvalue_diffexp.tsv"], contrast=config["diffexp"]["contrasts"]),
        expand(["results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt", "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"], contrast = config["diffexp"]["contrasts"], FC=config['FC'], adjp=config['adjp']),
        expand("results/diffexp/pairwise/{contrast}.diffexp.{adjp}.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"], adjp = config['adjp']),
        expand("results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf", contrast = config["diffexp"]["contrasts"]),
        expand(["results/diffexp/glimma-plots/{contrast}.ma_plot.html", "results/diffexp/glimma-plots/{contrast}.volcano_plot.html"],contrast = config["diffexp"]["contrasts"]),
        # "results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id),
        # expand(["results/diffexp/pairwise/enrichR/{contrast}-KEGG_2021_Human.upFC.{FC}.adjp.{adjp}.pdf", "results/diffexp/pairwise/enrichR/{contrast}-KEGG_2021_Human.downFC.{FC}.adjp.{adjp}.pdf"], contrast=config["diffexp"]["contrasts"], FC = config['FC'], adjp = config['adjp']),
        expand("results/diffexp/pairwise/enrichR/{contrast}.done", contrast=config["diffexp"]["contrasts"]),
        # expand("samples/star/{sample}_bam/{sample}_unmapped.fa", sample = SAMPLES),
        # expand("data/unmappedSeqs/{sample}_overRepseqCount.txt", sample = SAMPLES),
        "data/geneLengths.tsv",
        # expand("results/estSaturation/{plot}.png", plot = estSat_plotNames),
        # expand("results/estSaturation/{table}.tsv", table = estSat_tableNames),
        # expand("samples/bigwig/{sample}_cpm.bw", sample = SAMPLES),
        # expand("samples/bigwig/{sample}_fwd.bw", sample = SAMPLES),
        # expand("samples/bigwig/{sample}_rev.bw", sample = SAMPLES),
        "results/tables/compiled_readGC.tsv",


include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"
# include: "rules/unmapped.smk"
