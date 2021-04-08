rule collect_multiple_metrics:
    input:
         bam="results/alignment/recal/{sample}.bqsr.bam",
         ref=config['common']['genome']
    output:
        multiext("qc/{sample}",
                 ".alignment_summary_metrics",      # CollectAlignmentSummaryMetrics
                 ".insert_size_metrics",            # CollectInsertSizeMetrics
                 ".insert_size_histogram.pdf",      # CollectInsertSizeMetrics
                 ".quality_distribution_metrics",   # QualityScoreDistribution
                 ".quality_distribution.pdf",       # QualityScoreDistribution
                 ".gc_bias.detail_metrics",         # CollectGcBiasMetrics
                 ".gc_bias.summary_metrics",        # CollectGcBiasMetrics
                 ".gc_bias.pdf",                    # CollectGcBiasMetrics
                 ".quality_yield_metrics"           # CollectQualityYieldMetrics
                 )
    resources:
        mem_gb=4
    log:
        "logs/picard/multiple_metrics/{sample}.log"
    params:
        # optional parameters
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE "
    wrapper:
        "0.73.0/bio/picard/collectmultiplemetrics"

rule collect_wgs_metrics:
    input:
        bam="results/alignment/recal/{sample}.bqsr.bam",
        ref=config['common']['genome'],
    output:
        "results/qc/{sample}.wgs_metrics"
    log:
        "logs/picard/wgs_metrics/{sample}.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CollectWgsMetrics "
        "I={input.bam} "
        "O={output} "
        "R={input.ref} 2> {log}"

rule plot_wgs_insert:
    input:
        qcdir="results/qc",
        wgsfiles=expand("results/qc/{sample}.wgs_metrics", sample=samples.index),
        insertfiles=expand("results/qc/{sample}.insert_size_metrics", sample=samples.index)
    output:
        report("results/plots/qc/wgs_insert_metrics.pdf", caption="../report/qc.rst", category='QC'),
    params:
        samples=",".join(expand("{sample}", sample=samples.index)),
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript workflow/scripts/plot_picard_metrics.R "
        "--qcdir {input.qcdir} "
        "--samples {params.samples} "
        "--output {output}"
