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

configfile:"omic_config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz")
#SAMPLES = ["A10_ScreenNegative", "A1_Case", "A2_Case", "A3_Control", "A4_Control", "A5_Healthy", "A6_Healthy", "A7_InSitu", "A8_InSitu", "A9_ScreenNegative"]
#SAMPLES = [ "B10_ScreenNegative", "B1_Case", "B2_Case", "B3_Control", "B4_Control", "B5_Healthy", "B6_Healthy", "B7_InSitu", "B8_InSitu", "B9_ScreenNegative"]
#SAMPLES = ["C10_Control", "C1_ScreenNegative", "C2_ScreenNegative", "C3_InSitu", "C4_InSitu", "C5_Control", "C6_ScreenNegative", "C7_Case", "C8_Case", "C9_Healthy"]
#SAMPLES = ["D10_ScreenNegative", "D1_ScreenNegative", "D2_Case", "D3_InSitu", "D4_Healthy", "D5_Control", "D6_Control", "D7_Control", "D8_Case", "D9_Case"]
#SAMPLES = ["E10_Case", "E1_NA", "E2_InSitu", "E3_InSitu", "E4_ScreenNegative", "E5_ScreenNegative", "E6_Healthy", "E7_Healthy", "E8_Control", "E9_Control"]
#SAMPLES = ["F10_Case", "F1_Control", "F2_Control", "F3_NA", "F4_InSitu", "F5_InSitu", "F6_ScreenNegative", "F7_ScreenNegative", "F8_Healthy", "F9_Case"]
#SAMPLES = ["G10_Case", "G1_Control", "G2_Control", "G3_Case", "G4_Case", "G5_Healthy", "G6_InSitu", "G7_InSitu", "G8_NA", "G9_ScreenNegative"]
#SAMPLES = ["G4_Case", "G9_ScreenNegative", "H10_ScreenNegative", "H1_Control", "H2_ScreenNegative", "H3_Case", "H4_Control", "H5_Healthy", "H6_InSitu", "H7_Case", "H8_Control", "H9_Healthy", "I10_ScreenNegative", "I1_Case", "I2_InSitu", "I3_Case", "I4_Case", "I5_Healthy", "I6_InSitu", "I7_InSitu", "I8_Control", "I9_ScreenNegative"]
#SAMPLES = ["J10_Case", "J1_Control", "J2_Case", "J3_Case", "J4_Case", "J5_Control", "J6_ScreenNegative", "J7_ScreenNegative", "J8_Control", "J9_ScreenNegative", "K10_Case", "K1_Control", "K2_Case", "K3_Case", "K4_ScreenNegative", "K5_ScreenNegative", "K6_Case", "K7_Control", "K8_ScreenNegative", "K9_Control"]
#SAMPLES = ["L10_Control", "L1_ScreenNegative", "L2_Control", "L3_Case", "L4_Case", "L5_Case", "L6_ScreenNegative", "L7_ScreenNegative", "L8_InSitu", "L9_Control", "M10_Case", "M1_Control", "M2_ScreenNegative", "M3_ScreenNegative", "M4_ScreenNegative", "M5_Case", "M6_Healthy", "M7_Control", "M8_Case", "M9_Case"]
#SAMPLES = ["N10_Case", "N11_Case", "N1_Case", "N2_ScreenNegative", "N3_Healthy", "N4_Control", "N5_Case", "N6_Healthy", "N7_Case", "N8_ScreenNegative", "N9_Healthy"]


ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['R1','R2']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']
#plotNames = ['biotype_barplot', 'geneAttributes_barplot', 'totalReads_arranged']
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
        # expand("samples/htseq_count/{sample}_genecount.txt", sample = SAMPLES),
        # "data/{project_id}_genecounts.txt".format(project_id=config["project_id"]),
        # "data/{project_id}_genecounts_w_stats.txt".format(project_id=config["project_id"]),
        "data/{project_id}_genecounts.filt.txt".format(project_id=config["project_id"]),
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"]),
        # expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        # expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
        # expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        # expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        # "results/tables/read_coverage.txt",
        # expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        # expand("rseqc/geneBody_coverage/{sample}/{sample}.geneBodyCoverage.curves.pdf", sample = SAMPLES),
        # expand("results/readQC/{plot}.png", plot = readQC_plotNames),
        # expand("results/diffexp/pairwise/{contrast}.pca_plot.pdf", contrast = config["diffexp"]["contrasts"]),
        # "results/diffexp/group/LRT_pca.pdf",
        # "results/diffexp/group/MDS_table.txt",
        # "results/diffexp/group/LRT_density_plot.pdf",
        # expand(["results/diffexp/pairwise/{contrast}.qplot.pdf","results/diffexp/pairwise/{contrast}.qhist.pdf","results/diffexp/pairwise/{contrast}.qvalue_diffexp.tsv"],contrast=config["diffexp"]["contrasts"]),
        # expand(["results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt", "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt"], contrast = config["diffexp"]["contrasts"], FC=config['FC'], adjp=config['adjp']),
        # expand("results/diffexp/pairwise/{contrast}.diffexp.{adjp}.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"], adjp = config['adjp']),
        # expand("results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf", contrast = config["diffexp"]["contrasts"]),
        # expand(["results/diffexp/glimma-plots/{contrast}.ma_plot.html", "results/diffexp/glimma-plots/{contrast}.volcano_plot.html"],contrast = config["diffexp"]["contrasts"]),
        # "results/diffexp/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id),
        # expand("samples/star/{sample}_bam/{sample}_unmapped.fa", sample = SAMPLES),
        # expand("data/unmappedSeqs/{sample}_overRepseqCount.txt", sample = SAMPLES),
        # "data/geneLengths.tsv",
        # expand("results/estSaturation/{plot}.png", plot = estSat_plotNames),
        # expand("results/estSaturation/{table}.tsv", table = estSat_tableNames),
        # expand("samples/bigwig/{sample}_cpm.bw", sample = SAMPLES),
        # expand("samples/bigwig/{sample}_fwd.bw", sample = SAMPLES),
        # expand("samples/bigwig/{sample}_rev.bw", sample = SAMPLES),
        # "results/tables/compiled_readGC.tsv",


include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"
#include: "rules/unmapped.smk"
