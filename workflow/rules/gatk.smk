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
  wrapper:
    "0.73.0/bio/picard/markduplicates"

rule index_duplicates:
  input:
    "results/alignment/dedup/{sample}.bam"
  output:
    "results/alignment/dedup/{sample}.bam.bai"
  log:
    "logs/samtools/index/{sample}.log"
  wrapper:
    "0.73.0/bio/samtools/index"

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
    wrapper:
        "0.73.0/bio/gatk3/realignertargetcreator"

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
        java_temp=temp(directory("/tmp/gatk3_indelrealigner/{sample}")),
    log:
        "logs/gatk/indelrealigner/{sample}.log"
    params:
        extra=""  # optional
    threads: 8
    resources:
        mem_mb = 8192
    wrapper:
        "0.73.0/bio/gatk3/indelrealigner"

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
    wrapper:
        "0.73.0/bio/gatk3/baserecalibrator"

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
    wrapper:
        "0.73.0/bio/gatk3/printreads"
