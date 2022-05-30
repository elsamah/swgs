rule mark_duplicates:
  input:
    "results/alignment/mapped_reads/{sample}.sorted.bam"
  output:
    bam="results/alignment/dedup/{sample}.bam",
    metrics="results/alignment/dedup/{sample}.metric.txt"
  log:
    "logs/picard/dedup/{sample}.log"
  threads: 8
  params:
    "REMOVE_DUPLICATES=true",
    "ASSUME_SORT_ORDER='coordinate'"
  params:
    picard="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9b770440ff173434e53ee101c7452a0a/share/picard-2.26.0-0",
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    java -Xmx12g -jar {params.picard}/picard.jar \
    MarkDuplicates INPUT={input.bam} OUTPUT={output.dedup} \
    METRICS_FILE={output.metrics} \
    ASSUME_SORTED=true MAX_RECORDS_IN_RAM=100000 \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true
    """

rule index_duplicates:
  input:
    "results/alignment/dedup/{sample}.bam"
  output:
    "results/alignment/dedup/{sample}.bam.bai"
  log:
    "logs/samtools/index/{sample}.log"
  wrapper:
    "0.77.0/bio/samtools/index"

rule realigner_target_creator:
    input:
        bam="results/alignment/dedup/{sample}.bam",
        bai=rules.index_duplicates.output,
        ref=config['common']['genome'],
        known=get_snp_paths
    output:
        intervals="results/alignment/realign/{sample}.intervals",
        java_temp=temp(directory("gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk/indelrealigner/{sample}.realignertargetcreator.log",
    params:
        extra="", # optional
    resources:
        mem_mb=8192,
    threads: 8
    conda:
      "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
    shell:
    """
    gatk3 -Xmx8g -T RealignerTargetCreator \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -nt 4 \
    -I {input.bam} \
    -R {input.ref} \
    --interval_padding 100 \
    -known {input.known} \
    -dt None \
    -o {output.intervals}
    """

rule indelrealigner:
    input:
        bam="results/alignment/dedup/{sample}.bam",
        bai="results/alignment/dedup/{sample}.bam.bai",
        ref=config['common']['genome'],
        known=get_indel_paths,
        target_intervals="results/alignment/realign/{sample}.intervals"
    output:
        bam="results/alignment/realign/{sample}.bam",
        bai="results/alignment/realign/{sample}.bai",
#        java_temp=temp(directory("/tmp/gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk3/indelrealigner/{sample}.log"
    params:
        extra=""  # optional
    threads: 8
    conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
    resources:
        mem_mb = 8192
    shell:
    """
    gatk3 -Xmx12g -T IndelRealigner \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -I {input.bam} \
    -o {output} \
    -R {input.ref} \
    -targetIntervals {input.interval} \
    -known {input.known} \
    -dt None \
    -compress 0
    """

rule baserecalibrator:
    input:
        bam="results/alignment/realign/{sample}.bam",
        ref=config['common']['genome'],
        known=get_indel_paths
    output:
        "results/alignment/recal/{sample}.recal_data_table"
    log:
        "logs/gatk/bqsr/{sample}.recal.log",
    params:
        extra=combine_args(config["params"]["gatk"]["baserecalibrator"]),
    resources:
        mem_mb = 8192
    threads: 8
    conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
    shell:
    """
    gatk3 -Xmx12g -T BaseRecalibrator \
    -nct 4 \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -I {input.bam} \
    -o {output} \
    -R {input.ref} \
    -knownSites {input.known} \
    -rf BadCigar \
    -cov ReadGroupCovariate \
    -cov ContextCovariate \
    -cov CycleCovariate \
    -cov QualityScoreCovariate \
    -dt None
    """

rule printreads:
    input:
        bam="results/alignment/realign/{sample}.bam",
        ref=config['common']['genome'],
        recal_data="results/alignment/recal/{sample}.recal_data_table"
    output:
        "results/alignment/recal/{sample}.bqsr.bam"
    log:
        "logs/gatk/bqsr/{sample}.print.log"
    params:
        extra=combine_args(config["params"]["gatk"]["printreads"]),
    resources:
        mem_mb = 8192
    threads: 8
    conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
    wrapper:
    """
    gatk3 -Xmx12g -T PrintReads \
    {params.extra} \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -nct 4 \
    -I {input.bam} \
    -o {output} \
    -R {input.ref} \
    -BQSR {input.recal_data} \
    -rf BadCigar \
    -dt None
    """

rule symlink_bai:
    input:
        "results/alignment/recal/{sample}.bqsr.bam",
    params:
        "results/alignment/recal/{sample}.bqsr.bai",
    output:
        "results/alignment/recal/{sample}.bqsr.bam.bai",
    shell:
        "ln -s $(readlink -f {params}) {output}"
