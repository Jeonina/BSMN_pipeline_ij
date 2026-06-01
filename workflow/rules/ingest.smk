# =============================================================================
# Ingest rules: external alignment (BAM/CRAM) → calling-ready CRAM
#
# This stage replaces rules/mapping.smk when the input is a pre-aligned
# BAM/CRAM (samples.tsv bam schema). It normalizes any external alignment into
# the calling input anchor `results/mapping/{sample}/{sample}.cram` so
# rules/calling.smk and rules/filtering.smk run with ZERO changes.
#
# Validation (fail fast):
#   1. samtools quickcheck on the input.
#   2. Sequence-dictionary compatibility vs the reference .dict — a shared
#      contig with a different length aborts the run (scripts/validate_alignment_ref.py).
#
# Normalization:
#   * .cram input  → symlink to the anchor path + ensure .crai exists.
#   * .bam  input  → `samtools view -C -T REF` transcode to CRAM + index.
#
# All tools run inside the samtools Apptainer container (config/containers.yaml).
# =============================================================================

import os

# Absolute path to scripts/ (snakemake is always run from the project root).
_INGEST_SCRIPTS = os.path.abspath("scripts")

# Reference Picard dictionary used for the compatibility check.
_REF_DICT = config["ref"]["dict"]


# @MX:ANCHOR: [AUTO] ingest_alignment — produces the calling input anchor
#             results/mapping/{sample}/{sample}.cram consumed by mutect2_scatter,
#             get_pileup_summaries (calling.smk) and vaf_filter (filtering.smk).
# @MX:REASON: high fan_in contract — three downstream rules depend on this exact
#             output path + .crai sidecar. The reference-compatibility guard here
#             is the only barrier preventing a mis-referenced external alignment
#             from being called silently against the wrong genome.
rule ingest_alignment:
    """Normalize an external BAM/CRAM into the calling-ready CRAM anchor."""
    input:
        aln=get_external_alignment,
    output:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    params:
        ref=REF,
        ref_dict=_REF_DICT,
        samtools_sif=CONTAINERS["samtools"]["sif"],
        validate_script=os.path.join(_INGEST_SCRIPTS, "validate_alignment_ref.py"),
    log:
        "logs/ingest/{sample}/ingest_alignment.log",
    threads: max(2, workflow.cores // 4)
    resources:
        mem_mb=4000,
        runtime=480,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[ingest_alignment] START $(date -Iseconds)"
        echo "[ingest_alignment] sample={wildcards.sample}"
        echo "[ingest_alignment] input.aln={input.aln}"
        echo "[ingest_alignment] input.size=$(stat -c%s {input.aln} 2>/dev/null || echo unknown) bytes"
        echo "[ingest_alignment] ref={params.ref}"
        echo "[ingest_alignment] ref_dict={params.ref_dict}"
        echo "================================================================"

        mkdir -p $(dirname {output.cram})

        # 1. Structural sanity check on the input alignment.
        apptainer exec {params.samtools_sif} samtools quickcheck {input.aln}

        # 2. Sequence-dictionary compatibility vs the reference .dict.
        _header=$(dirname {output.cram})/.{wildcards.sample}.header.sam
        apptainer exec {params.samtools_sif} samtools view -H {input.aln} > "$_header"
        python {params.validate_script} --header "$_header" --dict {params.ref_dict}
        rm -f "$_header"

        # 3. Normalize to the calling input anchor (CRAM + .crai).
        case "{input.aln}" in
            *.cram)
                echo "[ingest_alignment] input is CRAM — linking to anchor"
                ln -sf "$(readlink -f {input.aln})" {output.cram}
                ;;
            *.bam)
                echo "[ingest_alignment] input is BAM — transcoding to CRAM"
                apptainer exec {params.samtools_sif} samtools view \
                    -@ {threads} -C -T {params.ref} \
                    -o {output.cram} {input.aln}
                ;;
            *)
                echo "[ingest_alignment] ERROR: unrecognized alignment extension" >&2
                exit 1
                ;;
        esac

        apptainer exec {params.samtools_sif} samtools index {output.cram}
        _exit=$?

        echo "================================================================"
        echo "[ingest_alignment] exit_code=$_exit"
        echo "[ingest_alignment] output.cram.size=$(stat -c%s {output.cram} 2>/dev/null || echo 0) bytes"
        echo "[ingest_alignment] END $(date -Iseconds)"
        echo "================================================================"
        """


rule samtools_flagstat_ingest:
    """QC: samtools flagstat on the ingested CRAM (matches mapping.smk target)."""
    input:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        flagstat="results/mapping/{sample}/flagstat.txt",
    params:
        ref=REF,
        samtools_sif=CONTAINERS["samtools"]["sif"],
    log:
        "logs/ingest/{sample}/flagstat.log",
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
