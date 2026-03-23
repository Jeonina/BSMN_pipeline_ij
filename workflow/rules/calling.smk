# =============================================================================
# Calling rules: CRAM → filtered somatic VCF
#
# Tumor-only Mutect2 with chromosome-level scatter-gather.
# All tools run inside Apptainer containers (config/containers.yaml).
# Chromosomes defined in config/config.yaml  calling.chromosomes.
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
        germline=_calling.get("germline_resource", ""),
        pon_flag=lambda wildcards: (
            f"--panel-of-normals {_calling.get('pon', {}).get('vcf', '')}"
            if _calling.get("pon", {}).get("use", False)
            else ""
        ),
        extra=_calling.get("mutect2_extra", ""),
        bqsr_mem=RESOLVED["bqsr_memory_gb"],
        tmpdir="results/calling/{sample}/tmp/{chrom}",
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/mutect2_scatter.{chrom}.log",
    threads: 2
    resources:
        mem_mb=lambda wildcards: RESOLVED["bqsr_memory_gb"] * 1024,
        runtime=2880,
    shell:
        """
        mkdir -p $(dirname {output.vcf}) {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx{params.bqsr_mem}G -Djava.io.tmpdir={params.tmpdir}" \
            Mutect2 \
            -R {params.ref} \
            -I {input.cram} \
            --tumor-sample {wildcards.sample} \
            --germline-resource {params.germline} \
            {params.pon_flag} \
            -L {wildcards.chrom} \
            --f1r2-tar-gz {output.f1r2} \
            -O {output.vcf} \
            --stats {output.stats} \
            {params.extra} 2> {log}
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
        tmpdir="results/calling/{sample}/tmp/merge",
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/merge_vcfs.log",
    threads: 1
    resources:
        mem_mb=4096,
        runtime=240,
    shell:
        """
        mkdir -p {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx4G -Djava.io.tmpdir={params.tmpdir}" \
            MergeVcfs \
            {params.vcf_flags} \
            -O {output.vcf} 2> {log}
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
        mem_mb=2048,
        runtime=60,
    shell:
        """
        apptainer exec {params.gatk_sif} \
            gatk MergeMutectStats \
            {params.stats_flags} \
            -O {output.stats} 2> {log}
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
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/learn_read_orientation.log",
    threads: 1
    resources:
        mem_mb=4096,
        runtime=120,
    shell:
        """
        apptainer exec {params.gatk_sif} \
            gatk LearnReadOrientationModel \
            {params.f1r2_flags} \
            -O {output.model} 2> {log}
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
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/get_pileup_summaries.log",
    threads: 1
    resources:
        mem_mb=8192,
        runtime=480,
    shell:
        """
        apptainer exec {params.gatk_sif} \
            gatk GetPileupSummaries \
            -R {params.ref} \
            -I {input.cram} \
            -V {params.variants} \
            -L {params.variants} \
            -O {output.table} 2> {log}
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
        mem_mb=2048,
        runtime=60,
    shell:
        """
        apptainer exec {params.gatk_sif} \
            gatk CalculateContamination \
            -I {input.pileup} \
            --tumor-segmentation {output.segmentation} \
            -O {output.contamination} 2> {log}
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
        tmpdir="results/calling/{sample}/tmp/filter",
        gatk_sif=CONTAINERS["gatk"]["sif"],
    log:
        "logs/calling/{sample}/filter_mutect_calls.log",
    threads: 1
    resources:
        mem_mb=8192,
        runtime=240,
    shell:
        """
        mkdir -p {params.tmpdir}
        apptainer exec {params.gatk_sif} \
            gatk --java-options "-Xmx8G -Djava.io.tmpdir={params.tmpdir}" \
            FilterMutectCalls \
            -R {params.ref} \
            -V {input.vcf} \
            --stats {input.stats} \
            --ob-priors {input.orientation} \
            --contamination-table {input.contamination} \
            {params.extra} \
            -O {output.vcf} 2> {log}
        rm -rf {params.tmpdir}
        """
