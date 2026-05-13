# =============================================================================
# Mapping rules: FASTQ → analysis-ready CRAM
#
# All tools run inside Apptainer containers (config/containers.yaml).
# Parameters sourced from config/resolved_params.yaml (via auto_params.py).
# Based on BSMN pipeline aln_1~5 — aln_4 (IndelRealigner) excluded per
# GATK4 best practices.
#
# DAG:
#   validate_fastq_pair (per RG, pre-flight)
#       → bwa_mem_sort (per RG)
#           → merge_bams (per sample)
#               → mark_duplicates
#                   → base_recalibrator
#                       → apply_bqsr → CRAM
#                           → samtools_flagstat
# =============================================================================


_MAPPING_SCRIPTS = os.path.abspath("scripts")


rule validate_fastq_pair:
    """
    Pre-flight FASTQ pair validation (M-FIX-001 Bug 1).

    Catches R1/R2 mismatch, truncated gzip, and count mismatch in seconds
    rather than allowing bwa to run for 4+ hours before sambamba crashes.
    Gates bwa_mem_sort via the sentinel file `.fastq_pair.ok`.
    """
    input:
        unpack(get_fastqs),
    output:
        ok=touch("results/mapping/{sample}/validation/{sample}.{rg}.fastq_pair.ok"),
    params:
        script=os.path.join(_MAPPING_SCRIPTS, "validate_fastq_pair.py"),
        mode="--quick",
    log:
        "logs/mapping/{sample}/validate_fastq_pair.{rg}.log",
    threads: 2
    resources:
        mem_mb=2048,
        runtime=30,
    shell:
        """
        python {params.script} --r1 {input.fq1} --r2 {input.fq2} {params.mode} \
            > {log} 2>&1
        """


rule bwa_mem_sort:
    """aln_1: BWA-MEM → sambamba view → sambamba sort  (per readgroup)"""
    input:
        unpack(get_fastqs),
        ok="results/mapping/{sample}/validation/{sample}.{rg}.fastq_pair.ok",
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


rule base_recalibrator_scatter:
    """
    aln_5 step 1a: GATK4 BaseRecalibrator (per-chromosome scatter).
    Runs in parallel across chromosomes for significantly faster execution.
    """
    input:
        bam="results/mapping/{sample}/{sample}.markduped.bam",
        bai="results/mapping/{sample}/{sample}.markduped.bai",
    output:
        table=temp("results/mapping/{sample}/bqsr_scatter/{chrom}.recal.table"),
    params:
        ref=REF,
        dbsnp=config["known_sites"]["dbsnp"],
        mills=config["known_sites"]["mills"],
        indels=config["known_sites"]["indels"],
        tmpdir="results/mapping/{sample}/tmp/bqsr_{chrom}",
        bqsr_mem=RESOLVED.get("gatk_memory_gb", 8),
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/mapping/{sample}/base_recalibrator.{chrom}.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: RESOLVED.get("gatk_memory_gb", 8) * 1024 + 512,
        runtime=1440,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[base_recalibrator_scatter] START $(date -Iseconds)"
        echo "[base_recalibrator_scatter] sample={wildcards.sample} chrom={wildcards.chrom}"
        echo "[base_recalibrator_scatter] input.bam={input.bam}"
        echo "[base_recalibrator_scatter] java_heap={params.bqsr_mem}G"
        echo "================================================================"

        mkdir -p $(dirname {output.table}) {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.bqsr_mem}G -Djava.io.tmpdir={params.tmpdir}" \
            BaseRecalibrator \
            -R {params.ref} \
            --known-sites {params.dbsnp} \
            --known-sites {params.mills} \
            --known-sites {params.indels} \
            -I {input.bam} \
            -L {wildcards.chrom} \
            -O {output.table}

        echo "================================================================"
        echo "[base_recalibrator_scatter] output.size=$(stat -c%s {output.table} 2>/dev/null || echo 0) bytes"
        echo "[base_recalibrator_scatter] END $(date -Iseconds)"
        echo "================================================================"

        rm -rf {params.tmpdir}
        """


rule gather_bqsr_reports:
    """
    aln_5 step 1b: Merge per-chromosome BQSR recalibration tables.
    Uses GATK GatherBQSRReports to combine scatter results.
    """
    input:
        tables=get_scattered_recal_tables,
    output:
        table=temp("results/mapping/{sample}/recal_data.table"),
    params:
        table_flags=lambda wildcards, input: " ".join(f"-I {t}" for t in input.tables),
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/mapping/{sample}/gather_bqsr_reports.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: RESOLVED.get("gatk_memory_gb", 8) * 1024,
        runtime=60,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[gather_bqsr_reports] START $(date -Iseconds)"
        echo "[gather_bqsr_reports] sample={wildcards.sample}"
        echo "[gather_bqsr_reports] input_tables=$(echo {input.tables} | wc -w)"
        echo "================================================================"

        apptainer exec {params.gatk_sif} \
            gatk GatherBQSRReports \
            {params.table_flags} \
            -O {output.table}

        echo "================================================================"
        echo "[gather_bqsr_reports] output.size=$(stat -c%s {output.table} 2>/dev/null || echo 0) bytes"
        echo "[gather_bqsr_reports] END $(date -Iseconds)"
        echo "================================================================"
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
