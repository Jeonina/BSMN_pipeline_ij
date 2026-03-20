# =============================================================================
# Mapping rules: FASTQ → analysis-ready CRAM
# Based on BSMN pipeline aln_1 ~ aln_5 (skipping aln_4 IndelRealigner)
# =============================================================================


rule bwa_mem_sort:
    """aln_1: BWA-MEM → sambamba view → sambamba sort (per readgroup)"""
    input:
        unpack(get_fastqs),
    output:
        bam=temp("results/mapping/{sample}/align/{sample}.{rg}.sorted.bam"),
        bai=temp("results/mapping/{sample}/align/{sample}.{rg}.sorted.bam.bai"),
    params:
        ref=REF,
        rg=r"@RG\tID:{sample}.{rg}\tSM:{sample}\tPL:illumina\tLB:{sample}\tPU:{rg}",
        bwa_threads=config["mapping"]["bwa_threads"],
        sort_threads=config["mapping"]["sort_threads"],
        sort_memory=config["mapping"]["sort_memory"],
    log:
        "logs/mapping/{sample}/bwa_mem_sort.{rg}.log",
    conda:
        "../envs/mapping.yaml"
    threads: config["mapping"]["bwa_threads"]
    resources:
        mem_mb=3000,
        runtime=1440,
    shell:
        """
        (bwa mem -M -t {params.bwa_threads} \
            -R '{params.rg}' \
            {params.ref} '{input.fq1}' '{input.fq2}' \
        | sambamba view -S -f bam -l 0 /dev/stdin \
        | sambamba sort -m {params.sort_memory} -t {params.sort_threads} \
            -o {output.bam} /dev/stdin) 2> {log}
        """


rule merge_bams:
    """aln_2: sambamba merge (per sample). Single RG → rename."""
    input:
        bams=get_sorted_bams,
    output:
        bam=temp("results/mapping/{sample}/{sample}.merged.bam"),
        bai=temp("results/mapping/{sample}/{sample}.merged.bam.bai"),
    log:
        "logs/mapping/{sample}/merge_bams.log",
    conda:
        "../envs/mapping.yaml"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=720,
    run:
        bams = input.bams
        if len(bams) == 1:
            shell("mv {bams[0]} {output.bam} && samtools index {output.bam} 2> {log}")
        else:
            shell("sambamba merge -t {threads} {output.bam} {bams} 2> {log}")


rule mark_duplicates:
    """aln_3: Picard MarkDuplicates"""
    input:
        bam="results/mapping/{sample}/{sample}.merged.bam",
    output:
        bam=temp("results/mapping/{sample}/{sample}.markduped.bam"),
        bai=temp("results/mapping/{sample}/{sample}.markduped.bai"),
        metrics="results/mapping/{sample}/markduplicates_metrics.txt",
    params:
        java_mem=config["mapping"]["markdup_memory"],
        tmpdir="results/mapping/{sample}/tmp",
    log:
        "logs/mapping/{sample}/mark_duplicates.log",
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        mem_mb=3000,
        runtime=1440,
    shell:
        """
        mkdir -p {params.tmpdir}
        picard -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmpdir} \
            MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -METRICS_FILE {output.metrics} \
            -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            -CREATE_INDEX true \
            -TMP_DIR {params.tmpdir} 2> {log}
        rm -rf {params.tmpdir}
        """


rule base_recalibrator:
    """aln_5 step1: GATK4 BaseRecalibrator"""
    input:
        bam="results/mapping/{sample}/{sample}.markduped.bam",
        bai="results/mapping/{sample}/{sample}.markduped.bai",
    output:
        table=temp("results/mapping/{sample}/recal_data.table"),
    params:
        ref=REF,
        dbsnp=config["known_sites"]["dbsnp"],
        mills=config["known_sites"]["mills"],
        indels=config["known_sites"]["indels"],
        tmpdir="results/mapping/{sample}/tmp",
    log:
        "logs/mapping/{sample}/base_recalibrator.log",
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        mem_mb=3000,
        runtime=1440,
    shell:
        """
        mkdir -p {params.tmpdir}
        gatk --java-options "-Xmx2G -Djava.io.tmpdir={params.tmpdir}" \
            BaseRecalibrator \
            -R {params.ref} \
            --known-sites {params.dbsnp} \
            --known-sites {params.mills} \
            --known-sites {params.indels} \
            -I {input.bam} \
            -O {output.table} 2> {log}
        rm -rf {params.tmpdir}
        """


rule apply_bqsr:
    """aln_5 step2: GATK4 ApplyBQSR → CRAM output"""
    input:
        bam="results/mapping/{sample}/{sample}.markduped.bam",
        bai="results/mapping/{sample}/{sample}.markduped.bai",
        table="results/mapping/{sample}/recal_data.table",
    output:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    params:
        ref=REF,
        tmpdir="results/mapping/{sample}/tmp",
    log:
        "logs/mapping/{sample}/apply_bqsr.log",
    conda:
        "../envs/mapping.yaml"
    threads: 2
    resources:
        mem_mb=3000,
        runtime=1440,
    shell:
        """
        mkdir -p {params.tmpdir}
        gatk --java-options "-Xmx2G -Djava.io.tmpdir={params.tmpdir}" \
            ApplyBQSR \
            -R {params.ref} \
            --bqsr-recal-file {input.table} \
            -I {input.bam} \
            -O /dev/stdout \
        | samtools view -@ {threads} -C -T {params.ref} \
            -o {output.cram} 2> {log}
        samtools index {output.cram}
        rm -rf {params.tmpdir}
        """


rule samtools_flagstat:
    """aln_5 step3: QC stats"""
    input:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        flagstat="results/mapping/{sample}/flagstat.txt",
    log:
        "logs/mapping/{sample}/flagstat.log",
    conda:
        "../envs/mapping.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
    shell:
        """
        samtools flagstat -@ {threads} {input.cram} > {output.flagstat} 2> {log}
        """
