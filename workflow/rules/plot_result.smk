# rule for ploting result

rule plot_kraken:
    """
    script to plot result for output of kraken
    """
    input:
        report = OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_kraken.report"
    output:
        plot = OUTDIR + "/{species}/PoolSeq_check/plots/{species}_{sex}_kraken.png"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_result_kraken.R"

rule plot_blast:
    """
    script to plot result for output of blast
    """
    input:
        blast = OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_blast_output.tsv"
    output:
        plot = OUTDIR + "/{species}/PoolSeq_check/plots/{species}_{sex}_blast.png"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_result_blast.R"

