contrast = get_contrast

rule deseq2_init_cov:
    input:
        counts = "data/{project_id}_genecounts.filt.txt".format(project_id=config["project_id"])
    output:
        rds="results/diffexp/covariate/pairwise/{contrast}_cov_all.rds",
        rld_out = "results/diffexp/covariate/pairwise/{contrast}_cov_rlog_dds.rds"
    params:
        md=config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        contrast = get_contrast,
        threads = get_deseq2_threads,
        covariate = config["covariate"]
    conda:
        "../envs/permutation.yaml"
    shell:
        """Rscript scripts/deseq2-init_cov.R \
        --countsFile={input.counts} \
        --outDir=results/diffexp/covariate/pairwise \
        --metaFile={params.md} \
        --sampleID={params.sample_id} \
        --linear_model={params.linear_model} \
        --contrast='{params.contrast}' \
        --threads={params.threads} \
        --covariate={params.covariate}"""


rule deseq2_pairwise_cov:
    input:
        rds="results/diffexp/covariate/pairwise/{contrast}_cov_all.rds",
        rld = "results/diffexp/covariate/pairwise/{contrast}_cov_rlog_dds.rds"
    output:
        table="results/diffexp/covariate/pairwise/{contrast}.diffexp.tsv",
        ma_plot="results/diffexp/covariate/pairwise/{contrast}.ma_plot.pdf",
        p_hist="results/diffexp/covariate/pairwise/{contrast}.phist_plot.pdf",
        heatmap_plot = "results/diffexp/covariate/pairwise/{contrast}.heatmap_plot.pdf",
        panel_ma = "results/diffexp/covariate/pairwise/{contrast}.panel_ma.pdf",
        var_heat = "results/diffexp/covariate/pairwise/{contrast}.variance_heatmap.pdf",
        pca_plot = "results/diffexp/covariate/pairwise/{contrast}.pca_plot.pdf"
    params:
        contrast=get_contrast,
        linear_model = config["linear_model"],
        pca_labels = config["pca"]["labels"],
        sample_id = config["linear_model"]
    conda:
        "../envs/deseq2.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2_pairwise_cov.R"


rule deseq2_group_cov:
    input:
        countsFile = "data/{project_id}_genecounts.filt.txt".format(project_id = project_id)
    output:
        pca="results/diffexp/covariate/group/LRT_pca.pdf",
        sd_mean_plot="results/diffexp/covariate/group/LRT_sd_mean_plot.pdf",
        distance_plot = "results/diffexp/covariate/group/LRT_distance_plot.pdf",
        heatmap_plot = "results/diffexp/covariate/group/LRT_heatmap_plot.pdf",
        rds="results/diffexp/covariate/group/LRT_all.rds",
        rld_out = "results/diffexp/covariate/group/LRT_rlog_dds.rds"
    params:
        pca_labels = config["pca"]["labels"],
        plotCols = format_plot_columns,
        samples = config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        LRT = config["diffexp"]["LRT"],
        colors = config['colors']['rcolorbrewer'],
        discrete = config['colors']['discrete'],
        covariate = config['covariate']
    conda:
        "../envs/deseq2.yaml"
    shell:
        """Rscript scripts/deseq2_group_cov.R \
        --countsFile={input.countsFile} \
        --outDir=results/diffexp/covariate/group \
        --metaFile={params.samples} \
        --sampleID={params.sample_id} \
        --linear_model={params.linear_model} \
        --plotCols={params.plotCols} \
        --LRT='{params.LRT}' \
        --pca_labels='{params.pca_labels}' \
        --colors='{params.colors}' \
        --discrete='{params.discrete}' \
        --covariate={params.covariate}
        """

rule deseq2_QC_cov:
    input:
        rld="results/diffexp/covariate/group/LRT_rlog_dds.rds",
        rds="results/diffexp/covariate/group/LRT_all.rds"
    output:
        mds_plot="results/diffexp/covariate/group/MDS_plot.pdf",
        mds_table="results/diffexp/covariate/group/MDS_table.txt",
        heatmap_plot="results/diffexp/covariate/group/Heatmap_all_genes.pdf",
        sd_plot="results/diffexp/covariate/group/stdev_plot.pdf",
        rlogCounts_plot="results/diffexp/covariate/group/rlog_counts_violinPlot.pdf",
        rlogCounts_fac_plot="results/diffexp/covariate/group/rlog_counts_faceted_violinPlot.pdf",
        counts_plot="results/diffexp/covariate/group/counts_violinPlot.pdf",
        counts_fac_plot="results/diffexp/covariate/group/counts_faceted_violinPlot.pdf"
    params:
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        plot_cols = format_plot_columns,
        colors = config['colors']['rcolorbrewer'],
        discrete = config['colors']['discrete']
    conda:
        "../envs/deseq2_QC.yaml"
    shell:
        """Rscript scripts/QC.R \
        --rld={input.rld} \
        --rds={input.rds} \
        --outDir=results/diffexp/covariate/group \
        --sampleID={params.sample_id} \
        --Type={params.linear_model} \
        --plot_cols={params.plot_cols} \
        --colors='{params.colors}' \
        --discrete='{params.discrete}'
        """

rule deseq2_qplot_cov:
    input:
        stats_table="results/diffexp/covariate/pairwise/{contrast}.diffexp.tsv",
    output:
        qplot="results/diffexp/covariate/pairwise/{contrast}.qplot.pdf",
        qhist="results/diffexp/covariate/pairwise/{contrast}.qhist.pdf",
        table = "results/diffexp/covariate/pairwise/{contrast}.qvalue_diffexp.tsv"
    params:
        contrast=get_contrast,
    conda:
        "../envs/qplot_env.yaml"
    script:
        "../scripts/qplot.R"


rule deseq2_density_cov:
    input:
        rld = "results/diffexp/covariate/group/LRT_rlog_dds.rds"
    output:
        density="results/diffexp/covariate/group/LRT_density_plot.pdf",
    params:
        linear_model = config["linear_model"],
        project_id = config["project_id"],
        colors = config['colors']['rcolorbrewer'],
        discrete = config['colors']['discrete']
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/density_plot.R"


# rule GO:
#     input:
#         degFile="results/diffexp/pairwise/{contrast}.diffexp.tsv"
#     output:
#         "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.downFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"]),
#         "results/diffexp/pairwise/GOterms/{{contrast}}.diffexp.upFC.{FC}.adjp.{adjp}_BP_GO.txt".format(FC = config["FC"],adjp=config["adjp"])
#     params:
#         contrast = get_contrast,
#         assembly = config["assembly"],
#         printTree = config["printTree"],
#         FC = config["FC"],
#         adjp = config["adjp"]
#     conda:
#         "../envs/runGO.yaml"
#     script:
#         "../scripts/runGOforDESeq2.R"


# rule volcano_cov:
#     input:
#         degFile="results/diffexp/covariate/pairwise/{contrast}.diffexp.tsv"
#     output:
#         volcano_plot="results/diffexp/covariate/pairwise/{{contrast}}.diffexp.{adjp}.VolcanoPlot.pdf".format(adjp=config["adjp"])
#     params:
#         contrast = get_contrast,
#         FC = config["FC"],
#         adjp = config["adjp"]
#     conda:
#         "../envs/deseq2.yaml"
#     script:
#         "../scripts/RNAseq_makeVolcano.R"


# rule permutation:
#     input:
#         counts = "data/{project_id}_genecounts.filt.txt".format(project_id=config["project_id"])
#     output:
#         numGenes = "results/diffexp/covariate/pairwise/permutationTest/{contrast}.number.diff.genes.csv",
#         permList = "results/diffexp/covariate/pairwise/permutationTest/{contrast}.permutation.list.csv",
#         histogram = "results/diffexp/covariate/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf"
#     params:
#         contrast = get_contrast,
#         samples = config["omic_meta_data"],
#         sample_id = config["sample_id"],
#         linear_model = config["linear_model"]
#     conda:
#         "../envs/permutation.yaml"
#     script:
#         "../scripts/permutation_test.R"


rule run_glimma_cov:
    input:
        rds="results/diffexp/covariate/pairwise/{contrast}_cov_all.rds"
    output:
        ma_plot = "results/diffexp/covariate/glimma-plots/{contrast}.ma_plot.html",
        volcano_plot = "results/diffexp/covariate/glimma-plots/{contrast}.volcano_plot.html",
    params:
        contrast = get_contrast,
        condition = config["linear_model"],
    conda:
        "../envs/glimma_env.yaml"
    script:
        "../scripts/run_glimma_cov.R"


rule run_glimma_mds_cov:
    input:
        rds="results/diffexp/covariate/group/LRT_all.rds"
    output:
        mds_plot = "results/diffexp/covariate/glimma-plots/{project_id}.mds_plot.html".format(project_id=project_id),
    params:
        project_id = config["project_id"],
    conda:
        "../envs/glimma_env.yaml"
    script:
        "../scripts/run_glimma_mds_cov.R"


rule runEnrichR_cov:
    input:
        deg="results/diffexp/covariate/pairwise/{contrast}.qvalue_diffexp.tsv"
    output:
        # "results/diffexp/pairwise/enrichR/{{contrast}}-KEGG_2021_Human.upFC.{FC}.adjp.{adjp}.pdf".format(FC = config["FC"], adjp = config['adjp']),
        # "results/diffexp/pairwise/enrichR/{{contrast}}-KEGG_2021_Human.downFC.{FC}.adjp.{adjp}.pdf".format(FC = config["FC"], adjp = config['adjp'])
        touch("results/diffexp/covariate/pairwise/enrichR/{contrast}.done")
    params:
        metaFile = config['omic_meta_data'],
        annoFile = config['filter_anno'],
        sampleID = config['sample_id'],
        padj = config['adjp'],
        FC = config['FC']
    conda:
        "../envs/runEnrichR.yaml"
    shell:
        """Rscript scripts/runEnrichR.R \
        --degFile={input.deg} \
        --metaFile={params.metaFile} \
        --annoFile={params.annoFile} \
        --outDir=results/diffexp/covariate/pairwise/enrichR \
        --sampleID={params.sampleID} \
        --padj={params.padj} \
        --FC={params.FC}
        """
