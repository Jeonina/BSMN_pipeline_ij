# =============================================================================
# Mapping rules: FASTQ → analysis-ready CRAM
#
# All tools run inside Apptainer containers (config/containers.yaml).
# Parameters sourced from config/resolved_params.yaml (via auto_params.py).
# Based on BSMN pipeline aln_1~5 — aln_4 (IndelRealigner) excluded per
# GATK4 best practices.
#
# DAG:
#   bwa_mem_sort (per RG)
#       → merge_bams (per sample)
#           → mark_duplicates
#               → base_recalibrator
#                   → apply_bqsr → CRAM
#                       → samtools_flagstat
# =============================================================================


rule bwa_mem_sort:
    """aln_1: BWA-MEM → sambamba view → sambamba sort  (per readgroup)"""
    input:
        unpack(get_fastqs),
    output:
        bam=temp("results/mapping/{sample}/align/{sample}.{rg}.sorted.bam"),
        bai=temp("results/mapping/{sample}/align/{sample}.{rg}.sorted.bam.bai"),
    params:
        ref=REF,
        rg=r"@RG\tID:{sample}.{rg}\tSM:{sample}\tPL:illumina\tLB:{sample}\tPU:{rg}",
        bwa_threads=lambda wildcards, threads: threads,
        sort_threads=lambda wildcards, threads: max(2, threads // 4),
        sort_memory=RESOLVED["sort_memory"],
        bwa_sif=CONTAINERS["bwa"]["sif"],
        sambamba_sif=CONTAINERS["sambamba"]["sif"],
    log:
        "logs/mapping/{sample}/bwa_mem_sort.{rg}.log",
    threads: max(4, workflow.cores // 2)
    resources:
        mem_mb=lambda wildcards, threads: threads * 3000,
        runtime=1440,
    shell:
        """
        mkdir -p $(dirname {output.bam})
        (apptainer exec {params.bwa_sif} bwa mem \
            -M -t {params.bwa_threads} \
            -R '{params.rg}' \
            {params.ref} '{input.fq1}' '{input.fq2}' \
        | apptainer exec {params.sambamba_sif} sambamba view \
            -S -f bam -l 0 /dev/stdin \
        | apptainer exec {params.sambamba_sif} sambamba sort \
            -m {params.sort_memory} -t {params.sort_threads} \
            -o {output.bam} /dev/stdin) 2> {log}
        """


rule merge_bams:
    """aln_2: sambamba merge (per sample).  Single RG → copy + index."""
    input:
        bams=get_sorted_bams,
    output:
        bam=temp("results/mapping/{sample}/{sample}.merged.bam"),
        bai=temp("results/mapping/{sample}/{sample}.merged.bam.bai"),
    params:
        sambamba_sif=CONTAINERS["sambamba"]["sif"],
        samtools_sif=CONTAINERS["samtools"]["sif"],
    log:
        "logs/mapping/{sample}/merge_bams.log",
    threads: max(2, workflow.cores // 8)
    resources:
        mem_mb=2000,
        runtime=720,
    run:
        bams = list(input.bams)
        log_f = str(log)
        if len(bams) == 1:
            # Single readgroup: copy instead of merge, then rebuild index
            shell(
                f"cp {bams[0]} {output.bam} 2>> {log_f}"
                f" && apptainer exec {params.samtools_sif}"
                f" samtools index {output.bam} 2>> {log_f}"
            )
        else:
            bam_list = " ".join(bams)
            # sambamba merge auto-creates the BAI
            shell(
                f"apptainer exec {params.sambamba_sif}"
                f" sambamba merge -t {threads} {output.bam} {bam_list} 2> {log_f}"
            )


rule mark_duplicates:
    """
    aln_3: Picard MarkDuplicates.
    OPTICAL_DUPLICATE_PIXEL_DISTANCE is auto-set per sequencer type
    (NovaSeq=2500, HiSeq=100) by scripts/auto_params.py.
    """
    input:
        bam="results/mapping/{sample}/{sample}.merged.bam",
    output:
        bam=temp("results/mapping/{sample}/{sample}.markduped.bam"),
        bai=temp("results/mapping/{sample}/{sample}.markduped.bai"),
        metrics="results/mapping/{sample}/markduplicates_metrics.txt",
    params:
        java_mem=RESOLVED["markdup_memory"],
        odpd=RESOLVED["optical_duplicate_pixel_distance"],
        tmpdir="results/mapping/{sample}/tmp",
        picard_sif=CONTAINERS["picard"]["sif"],
        picard_jar=CONTAINERS["picard"]["jar"],
    log:
        "logs/mapping/{sample}/mark_duplicates.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: (
            int(str(RESOLVED["markdup_memory"]).rstrip("G")) * 1024 + 512
        ),
        runtime=1440,
    shell:
        """
        mkdir -p {params.tmpdir}
        apptainer exec {params.picard_sif} \
            java -Xmx{params.java_mem} -Djava.io.tmpdir={params.tmpdir} \
            -jar {params.picard_jar} MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -METRICS_FILE {output.metrics} \
            -OPTICAL_DUPLICATE_PIXEL_DISTANCE {params.odpd} \
            -CREATE_INDEX true \
            -TMP_DIR {params.tmpdir} 2> {log}
        rm -rf {params.tmpdir}
        """


rule base_recalibrator:
    """aln_5 step 1: GATK4 BaseRecalibrator"""
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
        bqsr_mem=RESOLVED["bqsr_memory_gb"],
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/mapping/{sample}/base_recalibrator.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: RESOLVED["bqsr_memory_gb"] * 1024 + 512,
        runtime=1440,
    shell:
        """
        mkdir -p {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.bqsr_mem}G -Djava.io.tmpdir={params.tmpdir}" \
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
    """
    aln_5 step 2: GATK4 ApplyBQSR → CRAM via samtools.

    GATK writes recalibrated BAM to /dev/stdout; samtools converts to CRAM
    in the same pipe.  Both stderr streams are captured to the log file.
    """
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
        bqsr_mem=RESOLVED["bqsr_memory_gb"],
        gatk_sif=CONTAINERS["gatk"]["sif"],
        samtools_sif=CONTAINERS["samtools"]["sif"],
    log:
        "logs/mapping/{sample}/apply_bqsr.log",
    threads: max(2, workflow.cores // 8)
    resources:
        mem_mb=lambda wildcards: RESOLVED["bqsr_memory_gb"] * 1024 + 1024,
        runtime=1440,
    shell:
        """
        mkdir -p {params.tmpdir}
        (apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.bqsr_mem}G -Djava.io.tmpdir={params.tmpdir}" \
            ApplyBQSR \
            -R {params.ref} \
            --bqsr-recal-file {input.table} \
            -I {input.bam} \
            -O /dev/stdout \
        | apptainer exec {params.samtools_sif} samtools view \
            -@ {threads} -C -T {params.ref} \
            --output-fmt-option version=3.0 \
            -o {output.cram} -) 2> {log}
        apptainer exec {params.samtools_sif} \
            samtools index {output.cram} 2>> {log}
        rm -rf {params.tmpdir}
        """


rule samtools_flagstat:
    """QC: samtools flagstat on final CRAM"""
    input:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        flagstat="results/mapping/{sample}/flagstat.txt",
    params:
        ref=REF,
        samtools_sif=CONTAINERS["samtools"]["sif"],
    log:
        "logs/mapping/{sample}/flagstat.log",
    threads: min(4, max(1, workflow.cores // 16))
    resources:
        mem_mb=1000,
        runtime=60,
    shell:
        """
        apptainer exec {params.samtools_sif} samtools flagstat \
            -@ {threads} \
            --input-fmt-option reference={params.ref} \
            {input.cram} > {output.flagstat} 2> {log}
        """
