rule bwa_mem:
  input:
    reads=get_fastqs
  output:
    "results/alignment/mapped_reads/{sample}.sorted.bam"
  log:
    "logs/bwa_mem/{sample}.log"
  threads: 8
  params:
    index=config['params']['bwa']['mem']['genome'],
    extra=get_rgid,
    sort="samtools",
    sort_order="coordinate"
  wrapper:
    "0.77.0/bio/bwa/mem"


rule samtools_index:
  input:
    "results/alignment/mapped_reads/{sample}.sorted.bam"
  output:
    "results/alignment/mapped_reads/{sample}.sorted.bam.bai"
  params:
    ""
  wrapper:
    "0.77.0/bio/samtools/index"
