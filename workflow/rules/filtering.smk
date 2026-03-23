# =============================================================================
# Filtering rules: Mutect2-filtered VCF → final somatic SNV list
#
# Four sequential filters per PIPELINE.md:
#   1. accessibility_filter  — 1KG strict mask (keeps only mask base == 'P')
#   2. germline_filter       — removes gnomAD AF > 0.001 known germline variants
#   3. vaf_filter            — binomial test p < 1e-6 AND alt_count >= 5
#   4. pon_mask_filter       — IUPAC FASTA-based panel-of-normals mask
#
# DAG (per sample):
#   results/calling/{sample}/{sample}.filtered.vcf.gz
#       → accessibility_filter  → {sample}.accessible.txt (temp)
#       → germline_filter       → {sample}.germline_filtered.txt (temp)
#       → vaf_filter            → {sample}.vaf_filtered.txt (temp)
#       → pon_mask_filter       → {sample}.final.txt  ← final output
# =============================================================================

import os

# Safe reference to filtering config (rules parsed even when stage != filtering)
_filtering = config.get("filtering", {})

# Absolute path to scripts/ directory (works regardless of --directory argument)
_SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), "..", "scripts")
)


rule accessibility_filter:
    """
    Accessibility filter using the 1KG strict mask FASTA.
    Extracts PASS SNVs from the Mutect2-filtered VCF, queries the
    1KG mask with samtools faidx, and keeps only variants where
    the mask base at the variant position is 'P' (accessible).
    Output: 4-column text (chrom  pos  ref  alt).
    """
    input:
        vcf="results/calling/{sample}/{sample}.filtered.vcf.gz",
        tbi="results/calling/{sample}/{sample}.filtered.vcf.gz.tbi",
    output:
        txt=temp("results/filtering/{sample}/{sample}.accessible.txt"),
    params:
        mask=_filtering.get("mask_1kg", ""),
        samtools_sif=CONTAINERS["samtools"]["sif"],
        bcftools_sif=CONTAINERS["bcftools"]["sif"],
    log:
        "logs/filtering/{sample}/accessibility_filter.log",
    threads: 1
    resources:
        mem_mb=8192,
        runtime=480,
    shell:
        """
        mkdir -p $(dirname {output.txt})
        apptainer exec {params.bcftools_sif} \
            bcftools view -H -f PASS -v snps {input.vcf} \
            | cut -f1,2,4,5 \
            | awk -v samtools_sif="{params.samtools_sif}" \
                  -v mask="{params.mask}" '{{
                cmd="apptainer exec " samtools_sif \
                    " samtools faidx " mask \
                    " " $1 ":" $2 "-" $2 " 2>/dev/null | tail -n1";
                cmd | getline mask_base;
                close(cmd);
                if (mask_base == "P") print $1 "\\t" $2 "\\t" $3 "\\t" $4
            }}' > {output.txt} 2> {log}
        """


rule germline_filter:
    """
    Germline filter: remove known germline variants (gnomAD AF > 0.001).
    Reads the gnomAD SNP list and removes any variant whose
    (chrom:pos:ref:alt) key appears in the set.
    """
    input:
        txt="results/filtering/{sample}/{sample}.accessible.txt",
    output:
        txt=temp("results/filtering/{sample}/{sample}.germline_filtered.txt"),
    params:
        gnomad_snps=_filtering.get("gnomad", {}).get("snps", ""),
        script=os.path.join(_SCRIPTS, "germline_filter.py"),
    log:
        "logs/filtering/{sample}/germline_filter.log",
    threads: 1
    resources:
        mem_mb=4096,
        runtime=120,
    shell:
        """
        python {params.script} \
            --variants {params.gnomad_snps} \
            {input.txt} > {output.txt} 2> {log}
        """


rule vaf_filter:
    """
    VAF filter: keep somatic candidates passing binomial test + min alt count.
    Criterion: binom_test(alt_n, depth, p=0.5, alternative='less') < 1e-6
               AND alt_n >= 5
    Pileup performed directly on the sample CRAM via samtools mpileup
    (mapping quality >= 20, base quality >= 20).
    """
    input:
        txt="results/filtering/{sample}/{sample}.germline_filtered.txt",
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        txt=temp("results/filtering/{sample}/{sample}.vaf_filtered.txt"),
    params:
        ref=REF,
        p_threshold=_filtering.get("vaf", {}).get("p_binom_threshold", 1.0e-6),
        min_alt=_filtering.get("vaf", {}).get("min_alt_count", 5),
        samtools_sif=CONTAINERS["samtools"]["sif"],
        script=os.path.join(_SCRIPTS, "vaf_filter.py"),
    log:
        "logs/filtering/{sample}/vaf_filter.log",
    threads: 4
    resources:
        mem_mb=8192,
        runtime=480,
    shell:
        """
        python {params.script} \
            --cram {input.cram} \
            --ref {params.ref} \
            --samtools-sif {params.samtools_sif} \
            --min-mapq 20 \
            --min-baseq 20 \
            --p-threshold {params.p_threshold} \
            --min-alt {params.min_alt} \
            --threads {threads} \
            {input.txt} > {output.txt} 2> {log}
        """


rule pon_mask_filter:
    """
    PON mask filter using an IUPAC-encoded panel-of-normals FASTA.
    For each variant, queries the PON FASTA with samtools faidx and
    applies the IUPAC degenerate base lookup table:
      * → Pass (position not in PON)
      N → Fail (always masked)
      A/C/G/T → Fail if alt matches
      R/Y/S/W/K/M/B/D/H/V → Fail if alt is in the degenerate base set
    Produces the final, analysis-ready somatic SNV list.
    """
    input:
        txt="results/filtering/{sample}/{sample}.vaf_filtered.txt",
    output:
        txt="results/filtering/{sample}/{sample}.final.txt",
    params:
        pon_fasta=_filtering.get("pon_mask", {}).get("fasta", ""),
        samtools_sif=CONTAINERS["samtools"]["sif"],
        script=os.path.join(_SCRIPTS, "pon_mask_filter.py"),
    log:
        "logs/filtering/{sample}/pon_mask_filter.log",
    threads: 1
    resources:
        mem_mb=4096,
        runtime=120,
    shell:
        """
        python {params.script} \
            --pon-fasta {params.pon_fasta} \
            --samtools-sif {params.samtools_sif} \
            {input.txt} > {output.txt} 2> {log}
        """
