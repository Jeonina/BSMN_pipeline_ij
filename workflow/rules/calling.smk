# =============================================================================
# Calling rules: CRAM → filtered somatic VCF
#
# Tumor-only Mutect2 with chromosome-level scatter-gather.
# All tools run inside Apptainer containers (config/containers.yaml).
# Chromosomes defined in config/config.yaml  calling.chromosomes.
#
# Resources are dynamically allocated based on system capabilities
# (via config/resolved_params.yaml from scripts/auto_params.py).
#
# DAG (per sample):
#   mutect2_scatter × N_chroms
#       → merge_vcfs
#       → filter_mutect_calls  ← final output
#   mutect2_scatter × N_chroms
#       → merge_mutect_stats ──┐
#   mutect2_scatter × N_chroms  │
#       → learn_read_orientation ┤→ filter_mutect_calls
#   CRAM                        │
#       → get_pileup_summaries  │
#           → calculate_contamination ┘
# =============================================================================

# Safe reference to calling config — rules are parsed even when stage=mapping
_calling = config.get("calling", {})

# Dynamic resource helpers
_gatk_mem_gb = RESOLVED.get("gatk_memory_gb", 8)
_bqsr_mem_gb = RESOLVED.get("bqsr_memory_gb", 16)
_n_chroms = max(1, len(CHROMOSOMES))

# Cores per Mutect2 scatter job.  Mutect2's only parallel section is the
# PairHMM, and it stays single-threaded unless --native-pair-hmm-threads is
# passed.  Reserving cores/n_chroms without that flag (the previous behaviour)
# just idled the reservation, so the count is now explicit and forwarded to
# GATK.  PairHMM scaling flattens past ~4 threads; more cores are better spent
# on additional concurrent chromosomes.
_mutect2_threads = max(1, int(_auto(_calling.get("mutect2_threads"), "mutect2_threads", 4)))


rule mutect2_scatter:
    """
    Per-chromosome Mutect2 tumor-only call.
    Produces raw VCF, per-chrom stats, and f1r2 counts for orientation model.
    """
    input:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        vcf=temp("results/calling/{sample}/scatter/{chrom}.vcf.gz"),
        tbi=temp("results/calling/{sample}/scatter/{chrom}.vcf.gz.tbi"),
        stats=temp("results/calling/{sample}/scatter/{chrom}.vcf.gz.stats"),
        f1r2=temp("results/calling/{sample}/scatter/{chrom}.f1r2.tar.gz"),
    params:
        ref=REF,
        germline_flag=(
            f"--germline-resource {_calling.get('germline_resource', '')}"
            if _calling.get("germline_resource", "")
            else ""
        ),
        pon_flag=lambda wildcards: (
            f"--panel-of-normals {_calling.get('pon', {}).get('vcf', '')}"
            if _calling.get("pon", {}).get("use", False)
            else ""
        ),
        extra=_calling.get("mutect2_extra", ""),
        java_mem=_bqsr_mem_gb,
        tmpdir=scratch("{sample}", "mutect2_{chrom}"),
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/mutect2_scatter.{chrom}.log",
    threads: min(_mutect2_threads, workflow.cores)
    resources:
        mem_mb=lambda wildcards: _bqsr_mem_gb * 1024 + 1024,
        runtime=2880,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[mutect2_scatter] START $(date -Iseconds)"
        echo "[mutect2_scatter] sample={wildcards.sample} chrom={wildcards.chrom}"
        echo "[mutect2_scatter] input.cram={input.cram}"
        echo "[mutect2_scatter] input.cram.size=$(stat -c%s {input.cram} 2>/dev/null || echo unknown) bytes"
        echo "[mutect2_scatter] threads={threads}"
        echo "[mutect2_scatter] java_heap={params.java_mem}G"
        echo "[mutect2_scatter] mem_mb=$(({params.java_mem} * 1024))"
        echo "[mutect2_scatter] ref={params.ref}"
        echo "[mutect2_scatter] germline_flag={params.germline_flag}"
        echo "[mutect2_scatter] pon_flag={params.pon_flag}"
        echo "[mutect2_scatter] extra={params.extra}"
        echo "================================================================"

        mkdir -p $(dirname {output.vcf}) {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.java_mem}G -Djava.io.tmpdir={params.tmpdir}" \
            Mutect2 \
            -R {params.ref} \
            -I {input.cram} \
            --tumor-sample {wildcards.sample} \
            {params.germline_flag} \
            {params.pon_flag} \
            -L {wildcards.chrom} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2} \
            -O {output.vcf} \
            {params.extra}
        _exit=$?

        echo "================================================================"
        echo "[mutect2_scatter] exit_code=$_exit"
        echo "[mutect2_scatter] output.vcf.size=$(stat -c%s {output.vcf} 2>/dev/null || echo 0) bytes"
        echo "[mutect2_scatter] output.stats.size=$(stat -c%s {output.stats} 2>/dev/null || echo 0) bytes"
        echo "[mutect2_scatter] END $(date -Iseconds)"
        echo "================================================================"

        rm -rf {params.tmpdir}
        """


rule merge_vcfs:
    """Merge per-chromosome raw VCFs into a single file (picard MergeVcfs)."""
    input:
        vcfs=get_scattered_vcfs,
    output:
        vcf=temp("results/calling/{sample}/{sample}.merged.vcf.gz"),
        tbi=temp("results/calling/{sample}/{sample}.merged.vcf.gz.tbi"),
    params:
        vcf_flags=lambda wildcards, input: " ".join(f"-I {v}" for v in input.vcfs),
        java_mem=_gatk_mem_gb,
        tmpdir=scratch("{sample}", "merge"),
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/merge_vcfs.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024 + 512,
        runtime=240,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[merge_vcfs] START $(date -Iseconds)"
        echo "[merge_vcfs] sample={wildcards.sample}"
        echo "[merge_vcfs] input_vcf_count=$(echo {input.vcfs} | wc -w)"
        echo "[merge_vcfs] java_heap={params.java_mem}G"
        echo "[merge_vcfs] mem_mb=$(({params.java_mem} * 1024 + 512))"
        echo "================================================================"

        mkdir -p {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.java_mem}G -Djava.io.tmpdir={params.tmpdir}" \
            MergeVcfs \
            {params.vcf_flags} \
            -O {output.vcf}
        _exit=$?

        echo "================================================================"
        echo "[merge_vcfs] exit_code=$_exit"
        echo "[merge_vcfs] output.vcf.size=$(stat -c%s {output.vcf} 2>/dev/null || echo 0) bytes"
        echo "[merge_vcfs] END $(date -Iseconds)"
        echo "================================================================"

        rm -rf {params.tmpdir}
        """


rule merge_mutect_stats:
    """Merge per-chromosome Mutect2 stats (required by FilterMutectCalls)."""
    input:
        stats=get_scattered_stats,
    output:
        stats=temp("results/calling/{sample}/{sample}.merged.stats"),
    params:
        stats_flags=lambda wildcards, input: " ".join(f"--stats {s}" for s in input.stats),
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/merge_mutect_stats.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=60,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[merge_mutect_stats] START $(date -Iseconds)"
        echo "[merge_mutect_stats] sample={wildcards.sample}"
        echo "[merge_mutect_stats] input_stats_count=$(echo {input.stats} | wc -w)"
        echo "================================================================"

        apptainer exec {params.gatk_sif} \
            gatk MergeMutectStats \
            {params.stats_flags} \
            -O {output.stats}
        _exit=$?

        echo "================================================================"
        echo "[merge_mutect_stats] exit_code=$_exit"
        echo "[merge_mutect_stats] output.size=$(stat -c%s {output.stats} 2>/dev/null || echo 0) bytes"
        echo "[merge_mutect_stats] END $(date -Iseconds)"
        echo "================================================================"
        """


rule learn_read_orientation:
    """
    Build read-orientation artifact model from f1r2 counts.
    Passed to FilterMutectCalls as --ob-priors to correct FFPE/OxoG artifacts.
    """
    input:
        f1r2=get_scattered_f1r2,
    output:
        model=temp("results/calling/{sample}/{sample}.read_orientation_model.tar.gz"),
    params:
        f1r2_flags=lambda wildcards, input: " ".join(f"-I {f}" for f in input.f1r2),
        java_mem=_gatk_mem_gb,
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/learn_read_orientation.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=120,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[learn_read_orientation] START $(date -Iseconds)"
        echo "[learn_read_orientation] sample={wildcards.sample}"
        echo "[learn_read_orientation] input_f1r2_count=$(echo {input.f1r2} | wc -w)"
        echo "[learn_read_orientation] java_heap={params.java_mem}G"
        echo "================================================================"

        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.java_mem}G" \
            LearnReadOrientationModel \
            {params.f1r2_flags} \
            -O {output.model}
        _exit=$?

        echo "================================================================"
        echo "[learn_read_orientation] exit_code=$_exit"
        echo "[learn_read_orientation] output.size=$(stat -c%s {output.model} 2>/dev/null || echo 0) bytes"
        echo "[learn_read_orientation] END $(date -Iseconds)"
        echo "================================================================"
        """


rule get_pileup_summaries:
    """
    Pileup at common germline variant sites for contamination estimation.
    Uses the gnomAD small exac common variants resource.
    """
    input:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        table=temp("results/calling/{sample}/{sample}.pileup_summaries.table"),
    params:
        ref=REF,
        variants=_calling.get("contamination_resource", ""),
        intervals=" ".join(f"-L {c}" for c in CHROMOSOMES),
        java_mem=_gatk_mem_gb,
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/get_pileup_summaries.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024 + 512,
        runtime=480,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[get_pileup_summaries] START $(date -Iseconds)"
        echo "[get_pileup_summaries] sample={wildcards.sample}"
        echo "[get_pileup_summaries] input.cram={input.cram}"
        echo "[get_pileup_summaries] variants={params.variants}"
        echo "[get_pileup_summaries] java_heap={params.java_mem}G"
        echo "[get_pileup_summaries] chromosomes={params.intervals}"
        echo "================================================================"

        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.java_mem}G" \
            GetPileupSummaries \
            -R {params.ref} \
            -I {input.cram} \
            -V {params.variants} \
            {params.intervals} \
            -O {output.table}
        _exit=$?

        echo "================================================================"
        echo "[get_pileup_summaries] exit_code=$_exit"
        if [ -f {output.table} ]; then
            echo "[get_pileup_summaries] output_lines=$(wc -l < {output.table})"
            echo "[get_pileup_summaries] output.size=$(stat -c%s {output.table}) bytes"
        fi
        echo "[get_pileup_summaries] END $(date -Iseconds)"
        echo "================================================================"
        """


rule calculate_contamination:
    """Estimate sample contamination fraction from pileup summaries."""
    input:
        pileup="results/calling/{sample}/{sample}.pileup_summaries.table",
    output:
        contamination=temp("results/calling/{sample}/{sample}.contamination.table"),
        segmentation=temp("results/calling/{sample}/{sample}.tumor_segmentation.table"),
    params:
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/calculate_contamination.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=60,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[calculate_contamination] START $(date -Iseconds)"
        echo "[calculate_contamination] sample={wildcards.sample}"
        echo "[calculate_contamination] input.pileup={input.pileup}"
        echo "[calculate_contamination] input_lines=$(wc -l < {input.pileup})"
        echo "================================================================"

        apptainer exec {params.gatk_sif} \
            gatk CalculateContamination \
            -I {input.pileup} \
            --tumor-segmentation {output.segmentation} \
            -O {output.contamination}
        _exit=$?

        echo "================================================================"
        echo "[calculate_contamination] exit_code=$_exit"
        if [ -f {output.contamination} ]; then
            echo "[calculate_contamination] contamination_value=$(tail -1 {output.contamination} | cut -f2)"
        fi
        echo "[calculate_contamination] END $(date -Iseconds)"
        echo "================================================================"
        """


rule filter_mutect_calls:
    """
    Apply Mutect2 filters using orientation model + contamination estimates.
    Produces the final analysis-ready somatic VCF.
    """
    input:
        vcf="results/calling/{sample}/{sample}.merged.vcf.gz",
        tbi="results/calling/{sample}/{sample}.merged.vcf.gz.tbi",
        stats="results/calling/{sample}/{sample}.merged.stats",
        orientation="results/calling/{sample}/{sample}.read_orientation_model.tar.gz",
        contamination="results/calling/{sample}/{sample}.contamination.table",
    output:
        vcf="results/calling/{sample}/{sample}.filtered.vcf.gz",
        tbi="results/calling/{sample}/{sample}.filtered.vcf.gz.tbi",
    params:
        ref=REF,
        extra=_calling.get("filter_extra", ""),
        java_mem=_gatk_mem_gb,
        tmpdir=scratch("{sample}", "filter"),
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/filter_mutect_calls.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024 + 512,
        runtime=240,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[filter_mutect_calls] START $(date -Iseconds)"
        echo "[filter_mutect_calls] sample={wildcards.sample}"
        echo "[filter_mutect_calls] input.vcf={input.vcf}"
        echo "[filter_mutect_calls] input.stats={input.stats}"
        echo "[filter_mutect_calls] input.orientation={input.orientation}"
        echo "[filter_mutect_calls] input.contamination={input.contamination}"
        echo "[filter_mutect_calls] java_heap={params.java_mem}G"
        echo "[filter_mutect_calls] extra_flags={params.extra}"
        echo "================================================================"

        mkdir -p {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.java_mem}G -Djava.io.tmpdir={params.tmpdir}" \
            FilterMutectCalls \
            -R {params.ref} \
            -V {input.vcf} \
            --stats {input.stats} \
            --ob-priors {input.orientation} \
            --contamination-table {input.contamination} \
            {params.extra} \
            -O {output.vcf}
        _exit=$?

        echo "================================================================"
        echo "[filter_mutect_calls] exit_code=$_exit"
        echo "[filter_mutect_calls] output.vcf.size=$(stat -c%s {output.vcf} 2>/dev/null || echo 0) bytes"
        echo "[filter_mutect_calls] END $(date -Iseconds)"
        echo "================================================================"

        rm -rf {params.tmpdir}
        """
