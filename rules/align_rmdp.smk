rule trim_bbduk:
    input:
        fwd = "samples/raw/{sample}_R1.fastq.gz",
        rev = "samples/raw/{sample}_R2.fastq.gz"
    output:
        fwd = temp("samples/bbduk/{sample}/{sample}_R1_t.fastq.gz"),
        rev = temp("samples/bbduk/{sample}/{sample}_R2_t.fastq.gz"),
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
        temp("samples/bbduk/{sample}/good/{sample}_R1_t.good.fq.gz"),
        temp("samples/bbduk/{sample}/good/{sample}_R2_t.good.fq.gz"),
        temp("samples/bbduk/{sample}/bad/{sample}_R1_t.bad.fq.gz"),
        temp("samples/bbduk/{sample}/bad/{sample}_R2_t.bad.fq.gz"),
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
        "samples/star/{sample}_bam/Unmapped.out.mate1",
        "samples/star/{sample}_bam/Unmapped.out.mate2",
        "samples/star/{sample}_bam/Log.final.out"
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


rule filter_STARcounts:
    input:
        countsFile="data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    params:
        annoFile=config["filter_anno"],
        biotypes=config["biotypes"],
        mito=config['mito'],
        ercc=config['ERCC']
    shell:
        """Rscript scripts/RNAseq_filterCounts.R \
        --countsFile={input.countsFile} \
        --annoFile={params.annoFile} \
        --outDir=data \
        --biotypes={params.biotypes} \
        --mito={params.mito} \
        --ercc={params.ercc}"""


rule genecount:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/htseq_count/{sample}_genecount.txt"
    params:
        name = "genecount_{sample}",
        gtf = config["gtf_file"]
    conda:
        "../envs/omic_qc_wf.yaml"
    threads: 1
    shell:
        """htseq-count \
                -f bam \
                -r name \
                -s reverse \
                -m union \
                --additional-attr=gene_name \
                {input} \
                {params.gtf} > {output}"""


rule compile_counts:
    input:
        expand("samples/htseq_count/{sample}_genecount.txt", sample=SAMPLES)
    output:
        "data/{project_id}_genecounts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table.py"


rule compile_counts_and_stats:
    input:
        expand("samples/htseq_count/{sample}_genecount.txt",sample=SAMPLES)
    output:
        "data/{project_id}_genecounts_w_stats.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table_w_stats.py"


rule filter_genecounts:
    input:
        countsFile="data/{project_id}_genecounts.txt".format(project_id=config["project_id"])
    output:
        "data/{project_id}_genecounts.filt.txt".format(project_id=config["project_id"])
    params:
        annoFile=config["filter_anno"],
        biotypes=config["biotypes"],
        mito=config['mito'],
        ercc=config['ERCC']
    shell:
        """Rscript scripts/RNAseq_filterCounts.R \
        --countsFile={input.countsFile} \
        --annoFile={params.annoFile} \
        --outDir=data \
        --biotypes={params.biotypes} \
        --mito={params.mito} \
        --ercc={params.ercc}"""


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
        --outdir=results/readQC"""


rule cpm_tracks:
    input:
        bam = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        idx = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    output:
        wig = "samples/bigwig/{sample}_cpm.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM -bs 1 -b {input.bam} -o {output.wig}"


rule fwd_tracks:
    input:
        bam = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        idx = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    output:
       wig = "samples/bigwig/{sample}_fwd.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM --filterRNAstrand forward -bs 1 -b {input.bam} -o {output.wig}"


rule rev_tracks:
    input:
        bam = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        idx = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    output:
        wig = "samples/bigwig/{sample}_rev.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM --filterRNAstrand reverse -bs 1 -b {input.bam} -o {output.wig}"


rule make_geneLengthTable:
    output:
        "data/geneLengths.tsv"
    params:
        gtf = config["gtf_file"],
        anno = config["filter_anno"]
    conda:
        "../envs/estSaturation.yaml"
    shell:
        """Rscript scripts/make_geneLengthTable.R --gtfFile={params.gtf} --annoFile={params.anno} --outDir='data'"""


rule estSaturation:
    input:
        countsFile = "data/{project_id}_genecounts.txt".format(project_id=config["project_id"]),
        geneLengthsFile = "data/geneLengths.tsv"
    output:
        barplot = "results/estSaturation/barplot_nFeatures_bySample.png",
        facet1 = "results/estSaturation/facet_saturationCurve_byContrast.png",
        facet2 = "results/estSaturation/facet_saturationCurve_bySample.png",
        violin1 = "results/estSaturation/violin_nFeatureDistribution_byContrast.png",
        violin2 = "results/estSaturation/violin_saturationVariance_byContrast.png",
        satCurves = "results/estSaturation/saturationCurve_bySample.png",
        fpkm = "results/estSaturation/counts_fpkm.tsv",
        estSat = "results/estSaturation/estSaturation_results.tsv",
        summary = "results/estSaturation/summaryStatistics.tsv"
    params:
        md = config['omic_meta_data'],
        minCount = config['expression_threshold'],
        contrast = config['linear_model'],
        sampleID = config['sample_id']
    conda:
        "../envs/estSaturation.yaml"
    shell:
        """Rscript scripts/estSaturation.R \
        --countsFile={input.countsFile} \
        --geneLengthsFile={input.geneLengthsFile} \
        --mdFile={params.md} \
        --outDir=results/estSaturation \
        --minCount={params.minCount} \
        --contrast={params.contrast} \
        --sampleID={params.sampleID}"""
