rule trim_bbduk:
    input:
        fwd = "samples/raw/{sample}_R1.fastq.gz",
        rev = "samples/raw/{sample}_R2.fastq.gz"
    output:
        fwd = "samples/bbduk/{sample}/{sample}_R1_t.fastq.gz",
        rev = "samples/bbduk/{sample}/{sample}_R2_t.fastq.gz",
    params:
        ref=config["bb_adapter"]
    message:
        """--- Trimming."""
    conda:
        "../envs/bbmap.yaml"
    shell:
        """bbduk.sh -Xmx1g in1={input.fwd} in2={input.rev} out1={output.fwd} out2={output.rev} minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 ref={params.ref} hdist=1"""

rule afterqc_filter:
    input:
        fwd = "samples/bbduk/{sample}/{sample}_R1_t.fastq.gz",
        rev = "samples/bbduk/{sample}/{sample}_R2_t.fastq.gz"
    output:
        "samples/bbduk/{sample}/good/{sample}_R1_t.good.fq.gz",
        "samples/bbduk/{sample}/good/{sample}_R2_t.good.fq.gz",
        "samples/bbduk/{sample}/bad/{sample}_R1_t.bad.fq.gz",
        "samples/bbduk/{sample}/bad/{sample}_R2_t.bad.fq.gz",
        "samples/bbduk/{sample}/QC/{sample}_R1_t.fastq.gz.html",
        "samples/bbduk/{sample}/QC/{sample}_R1_t.fastq.gz.json",

    message:
        """---AfterQC"""
    conda:
        "../envs/afterqc.yaml"
    shell:
        """after.py -1 {input.fwd} -2 {input.rev} --report_output_folder=samples/bbduk/{wildcards.sample}/QC/ -g samples/bbduk/{wildcards.sample}/good/ -b samples/bbduk/{wildcards.sample}/bad/"""

rule fastqscreen:
    input:
        fwd = "samples/bbduk/{sample}/good/{sample}_R1_t.good.fq.gz",
        rev = "samples/bbduk/{sample}/good/{sample}_R2_t.good.fq.gz"
    output:
        "samples/fastqscreen/{sample}/{sample}_R1_t.good_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R1_t.good_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R1_t.good_screen.txt",
        "samples/fastqscreen/{sample}/{sample}_R2_t.good_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R2_t.good_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R2_t.good_screen.txt"
    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input.fwd} {input.rev}"""


rule fastqc:
    input:
        fwd = "samples/bbduk/{sample}/good/{sample}_R1_t.good.fq.gz",
        rev = "samples/bbduk/{sample}/good/{sample}_R2_t.good.fq.gz"
    output:
        fwd = "samples/fastqc/{sample}/{sample}_R1_t.good_fastqc.zip",
        rev = "samples/fastqc/{sample}/{sample}_R2_t.good_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """fastqc --outdir samples/fastqc/{wildcards.sample} --extract  -f fastq {input.fwd} {input.rev}"""

rule STAR:
    input:
        fwd = "samples/bbduk/{sample}/good/{sample}_R1_t.good.fq.gz",
        rev = "samples/bbduk/{sample}/good/{sample}_R2_t.good.fq.gz"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/{sample}_bam/Log.final.out",
        "samples/star/{sample}_bam/Unmapped.out.mate1",
        "samples/star/{sample}_bam/Unmapped.out.mate2"
    threads: 12
    params:
        gtf=config["gtf_file"]
    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix samples/star/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                --readFilesCommand zcat \
                --outReadsUnmapped Fastx \
                --twopassMode Basic
                """)

rule index:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/samtools_env.yaml"
    shell:
        """samtools index {input} {output}"""

rule star_statistics:
    input:
        expand("samples/star/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"


rule compile_star_counts:
    input:
        expand("samples/star/{sample}_bam/ReadsPerGene.out.tab",sample=SAMPLES)
    params:
        samples=SAMPLES
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_star_counts.py"


rule filter_counts:
    input:
        countsFile="data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    params:
        anno=config["filter_anno"],
        biotypes=config["biotypes"],
        mito=config['mito']
    script:
        "../scripts/RNAseq_filterCounts.R"


rule readQC:
    input:
        countsFile = "data/{project_id}_counts.txt".format(project_id=config["project_id"]),
        readDistFile = "results/tables/read_coverage.txt"
    output:
        readSummaryPlot = "results/readQC/totalReads_arranged.png",
        mappingPlot = "results/readQC/geneAttributes_barplot.png",
        biotypePlot = "results/readQC/biotype_barplot.png"
    params:
        annoFile = config['filter_anno'],
        metaFile = config['omic_meta_data'],
        contrast = config['linear_model'],
        corType = config['correlation_test'],
        plotCols = format_plot_columns
    conda:
        "../envs/readQC.yaml"
    shell:
        """Rscript scripts/RNAseq_readQC.R \
        --countsFile={input.countsFile} \
        --readDistFile={input.readDistFile} \
        --annoFile={params.annoFile} \
        --metaFile={params.metaFile} \
        --plotCols={params.plotCols} \
        --contrast={params.contrast} \
        --corType={params.corType} \
        --outdir=results/readQC_plots"""

