rule get_chromosomes:
    input:
        "results/alignment/recal/{sample}.bqsr.bam",
    output:
        "results/cnv/ichorcna/{sample}.chrs"
    params:
        read=1000,
        chrpattern="\"^chr[0-9XY]*$\"",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools idxstats {input} | "  # returns list of chr and reads
        "awk '$3 > {params.read}' - | " # remove chromosomes with <1000 reads
        "cut -f1 | "                    # selects chr column
        "grep {params.chrpattern} | "   # select chrs 1-22,X,Y
        "paste -s -d, - > {output}"     # comma-separated concatenation of chrs

rule readcounter:
    input:
        bam="results/alignment/recal/{sample}.bqsr.bam",
        bai="results/alignment/recal/{sample}.bqsr.bam.bai",
        chrs="results/cnv/ichorcna/{sample}.chrs",
    output:
        "results/cnv/ichorcna/{sample}.wig",
    params:
        window=config['params']['readcounter']['window'],
        quality=config['params']['readcounter']['quality'],
    log:
        "logs/cnv/ichorcna/{sample}_readcounter.log"
    conda:
        "../envs/r.yaml"
    shell:
        "readCounter "
        "--window {params.window} "
        "--quality {params.quality} "
        "--chromosome $(cat {input.chrs}) "
        "{input.bam} | "
        "sed \"s/chrom=chr/chrom=/\" "
        "> {output} 2> {log}"

rule get_RlibPath:
    output:
        "results/ref/libpath"
    conda:
        "../envs/r.yaml"
    priority:
        50
    shell:
        "Rscript -e \"cat(.libPaths(), '\n')\" > {output}"

rule ichorcna:
    input:
        bam="results/alignment/recal/{sample}.bqsr.bam",
        wig="results/cnv/ichorcna/{sample}.wig",
        rlib=rules.get_RlibPath.output,
        chrs=rules.get_chromosomes.output,
        ref=config['common']['genome'],
    output:
        dir=directory("results/cnv/ichorcna/{sample}"),
        file="results/cnv/ichorcna/{sample}/{sample}/{sample}_genomeWide.pdf",
    log:
        "logs/cnv/ichorcna/{sample}.log"
    conda:
        "../envs/r.yaml"
    params:
        ploidy=config['params']['ichorcna']['ploidy'],
        normal=config['params']['ichorcna']['normal'],
        maxCN=config['params']['ichorcna']['maxCN'],
        gc_wig=lambda w, input: get_ichorPath(input[2])['gc'],
        map_wig=lambda w, input: get_ichorPath(input[2])['map'],
        centromere=lambda w, input: get_ichorPath(input[2])['cen'],
        normal_panel=lambda w, input: get_ichorPath(input[2])['norm'],
        HOMD=config['params']['ichorcna']['include_HOMD'],
        chrs=lambda w, input: get_ichorChrs(input[3].format(sample=w.sample))['all'],
        chr_train=lambda w, input: get_ichorChrs(input[3].format(sample=w.sample))['train'],
        estimateNormal=config['params']['ichorcna']['estimateNormal'],
        estimatePloidy=config['params']['ichorcna']['estimatePloidy'],
        estimateScPrevalence=config['params']['ichorcna']['estimateScPrevalence'],
        sc_states=config['params']['ichorcna']['sc_states'],
        txnE=config['params']['ichorcna']['txnE'],
        txn_strength=config['params']['ichorcna']['txn_strength'],
        genome_style=config['params']['ichorcna']['genome_style'],
        genome_build=config['common']['build'],
    shell:
        "runIchorCNA.R "
        "--id {wildcards.sample} "
        "--WIG {input.wig} "
        "--ploidy '{params.ploidy}' "
        "--normal '{params.normal}' "
        "--maxCN {params.maxCN} "
        "--gcWig '{params.gc_wig}' "
        "--mapWig '{params.map_wig}' "
        "--centromere '{params.centromere}' "
        "--normalPanel '{params.normal_panel}' "
        "--includeHOMD {params.HOMD} "
        "--chrs '{params.chrs}' "
        "--chrTrain '{params.chr_train}' "
        "--estimateNormal {params.estimateNormal} "
        "--estimatePloidy {params.estimatePloidy} "
        "--estimateScPrevalence {params.estimateScPrevalence} "
        "--scStates '{params.sc_states}' "
        "--txnE {params.txnE} "
        "--txnStrength {params.txn_strength} "
        "--genomeStyle {params.genome_style} "
        "--outDir '{output.dir}' 2> {log}"

'''
        "--genomeBuild {params.genome_build} "
        "--NORMWIG {params.normalWig} "
        "--exons.bed {params.exons_bed} "
        "--minMapScore {params.min_map_score} "
        "--rmCentromereFlankLength {params.rmCentromereFlankLength} "
        "--coverage {params.coverage} "
        "--lambda {params.lambda} "
        "--lambdaScaleHyperParam {params.lambdaScaleHyperParam} "
        "--maxFracCNASubclone {params.maxFracCNASubclone} "
        "--maxFracGenomeSubclone {params.maxFracGenomeSubclone} "
        "--minSegmentBins {params.minSegmentBins} "
        "--altFracThreshold {params.altFracThreshold} "
        "--chrNormalize {params.chrNormalize} "
        "--genomeBuild {params.genomeBuild} "
        "--genomeStyle {params.genomeStyle} "
        "--normalizeMaleX {params.normalizeMalX} "
        "--fracReadsInChrYForMale {params.fracReadsInChrYForMale} "
        "--plotFileType {params.plotFileType} "
        "--plotYLim {params.plotYLim} "
        "--libdir {params.libdir} "
'''
