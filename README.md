![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.1-brightgreen.svg)[![Build Status](https://travis-ci.com/ohsu-cedar-comp-hub/Bulk-RNA-seq-pipeline-PE.svg?branch=master)](https://travis-ci.com/ohsu-cedar-comp-hub/Bulk-RNA-seq-pipeline-PE)
# Bulk-RNA-seq-pipeline-PE

Pipeline to run basic bulk RNA-seq analysis on paired-end data.

This is a package of Python and R scripts that enable reading, processing and analysis of -omics' datasets. 
This package implements the Snakemake management workflow system and is currently implemented to work with the cluster management and job scheduling system SLURM. 
This snakemake workflow utilizes conda installations to download and use packages for further analysis, so please ensure that you have installed miniconda prior to use.


Questions/issues
======================

Please add an issue to the repository. We would appreciate if your issue included sample code/files (as appropriate) so that we can reproduce your bug/issue.


Contributing
======================

We welcome contributors! For your pull requests, please include the following:

* Sample code/file that reproducibly causes the bug/issue
* Documented code providing fix
* Unit tests evaluating added/modified methods. 


Use
======================

Clone this repository into your working directory.

```
$ cd /path/to/working/directory
$ git clone https://github.com/ohsu-cedar-comp-hub/Bulk-RNA-seq-pipeline-PE.git
```

Go into the new directory -- this is now your working directory (wdir). 

```
$ cd Bulk-RNA-seq-pipeline-PE
```

Create a `samples/raw` directory.

```
$ mkdir -p samples/raw
```

Symbolically link the fastq files of your samples to the `wdir/samples/raw` directory using a bash script loop in your terminal. NOTE: you may need to change the file names to something meaningful to you.

```
ls -1 /path/to/data/archive/*.fastq.gz | while read fastq ; do ln -s $fastq samples/raw ; done
```

Upload your metadata file to the `data` directory, with the correct formatting:
* Columns should read:
```SampleID  Column2   Column3   ...```
* Each row should be a sample, with subsequent desired information provided in the columns (RNA extraction date, alias, time, etc.)
* File is tab-separated (.tsv)

Edit `omic_config.yaml` based on the project specifics.
* Project-specific file paths
    * `gtf_file`, `bed_file`, `star_index`, `filter_anno`, `omic_meta_data`
* Project details and specifications
    * `project_id`, `assembly`
* Parameters to select
    * `biotypes`, `mito`, `printTree`, `FC`, `adjp`, `seq_layout`, `correlation_test`, `expression_threshold`
* `DESeq2` details
    * `linear_model`, `sample_ID`, `LRT`
    * Base these off your metadata.tsv: `meta_columns_to_plot`, `pca` (labels)
    * Add appropriate contrasts based on your samples under `diffexp`
* Optional details
    * `colors`

Do a dry-run of snakemake to ensure proper execution before submitting it to the cluster (in your wdir).

```
$ conda activate <your_snakemake_environment>
$ snakemake -np --verbose
```

Once your files are symbolically linked, you can submit the job to exacloud via your terminal window.

```
$ sbatch submit_snakemake.sh
```

To see how the job is running, look at your queue.

```
$ squeue -u your_username
```

Detailed Workflow
=================================

Alignment
======================
1) Trimming & Filtering
    * Trimming and filtering of paired-end reads is performed using the trimming tool `bbduk`. For more information visit: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
    * Output files are located as `samples/bbduk/{sample}/{sample}_R*_t.good.fastq.gz`
2) Quality Analysis
    * Trimmed and filtered reads are assessed for quality using `AfterQC`, `FastQC` and `FastQScreen`. Output summaries are compiled into an interactive report by `MultiQC` at `./multiqc/{project_id}_QC.html`.
3) Alignment
    * Trimmed and filtered reads are aligned to the specified reference assembly using `STAR`
        * We included a two pass mode flag in order to increase the number of aligned reads
        * Unmapped reads are also saved within their respective directories
        * Output files are placed in `samples/star/{sample}_bam/` and include:
            * `Aligned.sortedByCoord.out.bam`
            * `ReadsPerGene.out.tab`
            * `Log.final.out`
    * We extracted the statistics from the `STAR` run and placed them in a table, summarizing the results across all samples from the `Log.final.out` output of STAR
        * Output is `results/tables/{project_id}_STAR_mapping_statistics.txt`
    * Any unmapped reads are saved as a FASTQ for each read pair. We also compiled these together into a single FASTA file.
        * R1 FASTQ: `samples/star/{sample}_bam/Unmapped.out.mate1`
        * R2 FASTQ: `samples/star/{sample}_bam/Unmapped.out.mate2`
        * Sample FASTA: `samples/star/{sample}_bam/{sample}_unmapped.fa`
4) Gene counts matrix
    * The gene counts are extracted and compiled into one matrix for each sample from the `STAR` counts
    * The raw, pre-filtered output is located as `data/{project_id}_counts.txt`
    * Based on the `omic_config.yaml` parameters `anno`, `biotypes`, and `mito`, the filtered gene counts matrix is located as `data/{project_id}_counts.filt.txt`

Additional Quality Analysis / Quality Check
======================
1) Read quality using `RSEQC`
    * `RSEQC` was used to check the quality of the reads by using a collection of commands from the `RSEQC` package:
        * Insertion Profile
        * Inner Distance
        * Clipping Profile
        * Read distribution
        * Read GC
        * Genebody coverage
    * For more information on these, visit: http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/rseqc/_build/html/index.html#usage-information
    * Output directory: `rseqc/`
2) Additional read quality metrics are stored in `results/readQC` and include:
    * Gene attribute fractions
    * RNA biotype distributions
    * Pairwise sample correlations
3) Sequencing saturation is estimated using `RNAseQC`. Results can be found in `results/estSaturation` and include: 
    * FPKM counts table & summary statistics
    * Sample feature quantifications
    * Sequencing saturation curves

Differential Expression Analysis (`DESeq2`)
======================
1) 

<!-- 2) QA/QC scripts to analyze the data as a whole 
    * The purpose of this analysis is to identify potential batch effects and outliers in the data
    * The outputs are located in `results/diffexp/` and are distributed in each analysis (i.e., `group` and `pairwise`)

    * The outputs to this are located in the `results` directory, and are distributed amongst 4 subdirectories, numbered `1 through 4`
        * `1`
            * A *boxplot* of the raw log2-transformed gene counts across all samples
            * A *boxplot* of the loess-transformed gene counts across all samples
            * A *scatter plot* comparing raw gene counts to loess-transformed gene counts
            * A *density plot* of raw log2-transformed gene counts across all samples 
            * A *density plot* of loess-transformed gene counts across all samples
            * A *scatter plot* of the standard deviation of raw log2-transformed gene counts across all samples
            * A *scatter plot* of the standard deviation of loess-transformed gene counts across all samples
        * `2`
            * A *heatmap* of all raw log2-transformed gene counts across samples
            * A *heatmap* of all loess-transformed gene counts across samples
                * These are generated to look for any batch effects in the data, due to date of extraction, or other factors
            * An *MDS Plot* for all samples, generated with the raw log2-transformed gene counts
            * An *MDS Plot* for all samples, generated with the loess-transformed gene counts
                * These are generated to look for outliers in the data
        * `3`
            * *p-value histograms* for each contrast specified in the `omic_config.yaml`
            * *q-value QC plot arrays* for each contrast specified in the `omic_config.yaml`
        * `4`
            * A *Heatmap* which looks at genes with a high FC and low q-value (very significant)
                * Takes genes with a FC>1.3, and ranks those by q-value. From this, a heatmap is generated for the top *50, 100 and 200* genes in this list
            * An *MDS Plot* which looks at the same subsets of genes as the Heatmap described above -->
            
Differential Expression Analysis (DESeq2)
======================
1) Initializing the DESeq2 object
    * Here, we run `DESeq2` on the genecounts table, which generates an RDS object and rlog counts
        * This includes the DE analysis across all samples
        * Output is located in the `results/diffexp/` 
    * From the dds object generated, we extract the normalized counts and generate a table with the results
        * Output is `results/tables/{project_id}_normed_counts.txt`
2) Generating plots
    * From the RDS object, we generate a collection of informative plots. These include:
        * *PCA Plot*
        * *Standard Deviation from the Mean Plot*
        * *Heatmap*
        * *Variance Heatmap*
        * *Distance Plot*
3) Differential Expression Analysis
    * We perform Differential Expression (DE) analysis for each contrast listed in the `omic_config.yaml`
    * Our output consists of DE gene count tables and a variety of plots
        * A table is generated for genes that are differentially expressed for each contrast
            * The output is placed in `results/diffexp/{contrast}.diffexp.tsv`
        * *MA Plots* are generated for each contrast
        * *p-histograms* are generated for each contrast
4) Differential Expression Plots
    * We use the output from DESeq2 to generate two types of plots:
        * Gene Ontology (GO) plots:
            * A `tree graph` describing the GO ID relationship for significantly up/downregulated genes in a given comparison
                * Output is located in `results/diffexp/GOterms`
            * A `bar graph` describing the enrichment and significance of GO IDs for up/downregulated genes in a given comparison
        * Volcano plots:
            * A `volcano plot` describing the distribution of up/downregulated genes in a given comparison
                * Output is located in `results/diffexp`
