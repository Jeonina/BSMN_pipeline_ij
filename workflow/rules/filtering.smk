# =============================================================================
# Filtering rules: Mutect2-filtered VCF → final somatic SNV list
#
# Five sequential filters per PIPELINE.md:
#   1. accessibility_filter  — 1KG strict mask (keeps only mask base == 'P')
#   2. germline_filter       — removes gnomAD AF > 0.001 known germline variants
#   3. vaf_filter            — binomial test p < 1e-6 AND alt_count >= 5
#   4. mayo_filter           — BSMN E-step: strand bias + repeat + multiallelic
#   5. pon_mask_filter       — IUPAC FASTA-based panel-of-normals mask
#
# Resources are dynamically allocated based on system capabilities
# (via config/resolved_params.yaml from scripts/auto_params.py).
#
# DAG (per sample):
#   results/calling/{sample}/{sample}.filtered.vcf.gz
#       → accessibility_filter  → {sample}.accessible.txt (temp)
#       → germline_filter       → {sample}.germline_filtered.txt (temp)
#       → vaf_filter            → {sample}.vaf_filtered.txt (temp)
#       → mayo_filter           → {sample}.mayo_filtered.txt (temp)
#       → pon_mask_filter       → {sample}.final.txt  ← final output
# =============================================================================

import os

# Safe reference to filtering config (rules parsed even when stage != filtering)
_filtering = config.get("filtering", {})

# Absolute path to scripts/ directory (snakemake is always run from project root)
_SCRIPTS = os.path.abspath("scripts")

# Dynamic resource helpers
_gatk_mem_gb = RESOLVED.get("gatk_memory_gb", 8)


rule accessibility_filter:
    """
    Accessibility filter using the 1KG strict mask BED.
    Extracts PASS SNVs from the Mutect2-filtered VCF and keeps only
    variants that fall within accessible regions in the BED file.
    Output: 4-column text (chrom  pos  ref  alt).
    """
    input:
        vcf="results/calling/{sample}/{sample}.filtered.vcf.gz",
        tbi="results/calling/{sample}/{sample}.filtered.vcf.gz.tbi",
    output:
        txt=temp("results/filtering/{sample}/{sample}.accessible.txt"),
    params:
        bed=_filtering.get("mask_1kg", ""),
        bcftools_sif=CONTAINERS["bcftools"]["sif"],
        script=os.path.join(_SCRIPTS, "accessibility_filter.py"),
    log:
        "logs/filtering/{sample}/accessibility_filter.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=480,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[accessibility_filter] START $(date -Iseconds)"
        echo "[accessibility_filter] sample={wildcards.sample}"
        echo "[accessibility_filter] input.vcf={input.vcf}"
        echo "[accessibility_filter] bed={params.bed}"
        echo "================================================================"

        mkdir -p $(dirname {output.txt})

        python {params.script} \
            --vcf {input.vcf} \
            --bed {params.bed} \
            --bcftools-sif {params.bcftools_sif} \
            > {output.txt}

        _kept=$(wc -l < {output.txt})
        echo "[accessibility_filter] variants_kept=$_kept"
        echo "[accessibility_filter] END $(date -Iseconds)"
        echo "================================================================"
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
        af_threshold=_filtering.get("gnomad", {}).get("af_threshold", 0.001),
        script=os.path.join(_SCRIPTS, "germline_filter.py"),
    log:
        "logs/filtering/{sample}/germline_filter.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
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
    Criterion: binom_test(alt_n, depth, p=0.5, alternative='less') < p_threshold
               AND alt_n >= min_alt
    Pileup performed directly on the sample CRAM via samtools mpileup.
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
        min_mapq=_filtering.get("vaf", {}).get("min_mapq", 20),
        min_baseq=_filtering.get("vaf", {}).get("min_baseq", 20),
        samtools_sif=CONTAINERS["samtools"]["sif"],
        script=os.path.join(_SCRIPTS, "vaf_filter.py"),
    log:
        "logs/filtering/{sample}/vaf_filter.log",
    threads: max(2, workflow.cores // 4)
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=480,
    shell:
        """
        python {params.script} \
            --cram {input.cram} \
            --ref {params.ref} \
            --samtools-sif {params.samtools_sif} \
            --min-mapq {params.min_mapq} \
            --min-baseq {params.min_baseq} \
            --p-threshold {params.p_threshold} \
            --min-alt {params.min_alt} \
            --threads {threads} \
            {input.txt} > {output.txt} 2> {log}
        """


rule mayo_filter:
    """
    Mayo filter — BSMN E.mayo_filters.sh step (strand bias + repeat + multiallelic).

    Runs the original strand_bias_filter.py and repeat_filter.py on the
    VAF-passing candidates, then mayo_filter.py keeps a candidate only when:
      repeat:         repeat_n < max_repeat_n AND repeat_length < max_repeat_length
      multiallelic:   total == ref_n + alt_n
      both strands:   alt_fwd >= 1 AND alt_rev >= 1
      strand balance: p_poisson >= min_strand_p OR p_fisher >= min_strand_p

    samtools (needed by both upstream scripts, which call it from PATH) is
    provided via a tiny PATH shim that delegates to the Apptainer container.
    The sample CRAM is read by injecting '-f {ref} {cram}' into samtools
    mpileup through strand_bias_filter.py's --bam argument.
    """
    input:
        txt="results/filtering/{sample}/{sample}.vaf_filtered.txt",
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        txt=temp("results/filtering/{sample}/{sample}.mayo_filtered.txt"),
    params:
        ref=REF,
        strand_tsv="results/filtering/{sample}/{sample}.strand.tsv",
        repeat_tsv="results/filtering/{sample}/{sample}.repeat.tsv",
        binshim="results/filtering/{sample}/_shim",
        max_repeat_n=_filtering.get("mayo", {}).get("max_repeat_n", 4),
        max_repeat_length=_filtering.get("mayo", {}).get("max_repeat_length", 10),
        min_strand_p=_filtering.get("mayo", {}).get("min_strand_p", 0.05),
        min_mapq=_filtering.get("mayo", {}).get("min_mapq", 20),
        min_baseq=_filtering.get("mayo", {}).get("min_baseq", 20),
        samtools_sif=CONTAINERS["samtools"]["sif"],
        strand_script=os.path.join(_SCRIPTS, "strand_bias_filter.py"),
        repeat_script=os.path.join(_SCRIPTS, "repeat_filter.py"),
        mayo_script=os.path.join(_SCRIPTS, "mayo_filter.py"),
    log:
        "logs/filtering/{sample}/mayo_filter.log",
    threads: max(2, workflow.cores // 4)
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=480,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[mayo_filter] START $(date -Iseconds)"
        echo "[mayo_filter] sample={wildcards.sample}"
        echo "[mayo_filter] input.txt={input.txt}"
        echo "================================================================"

        # samtools PATH shim → Apptainer container (both upstream scripts call
        # bare `samtools`; repeat_filter.py uses faidx, strand_bias_filter.py mpileup)
        SIF="{params.samtools_sif}"
        mkdir -p {params.binshim}
        cat > {params.binshim}/samtools <<SHIM
#!/bin/bash
exec apptainer exec $SIF samtools "\$@"
SHIM
        chmod +x {params.binshim}/samtools
        export PATH="$(cd {params.binshim} && pwd):$PATH"

        # strand bias (CRAM read via injected '-f {params.ref} {input.cram}')
        python {params.strand_script} \
            -b "-f {params.ref} {input.cram}" \
            -q {params.min_mapq} -Q {params.min_baseq} -n {threads} \
            {input.txt} > {params.strand_tsv}

        # repeat / STR context
        python {params.repeat_script} \
            -r {params.ref} -n {threads} \
            {input.txt} > {params.repeat_tsv}

        # join + apply the four mayo criteria
        python {params.mayo_script} \
            --strand {params.strand_tsv} \
            --repeat {params.repeat_tsv} \
            --max-repeat-n {params.max_repeat_n} \
            --max-repeat-length {params.max_repeat_length} \
            --min-strand-p {params.min_strand_p} \
            > {output.txt}

        rm -rf {params.binshim}
        _kept=$(wc -l < {output.txt})
        echo "[mayo_filter] variants_kept=$_kept"
        echo "[mayo_filter] END $(date -Iseconds)"
        echo "================================================================"
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
        txt="results/filtering/{sample}/{sample}.mayo_filtered.txt",
    output:
        txt="results/filtering/{sample}/{sample}.final.txt",
    params:
        pon_fasta=_filtering.get("pon_mask", {}).get("fasta", ""),
        samtools_sif=CONTAINERS["samtools"]["sif"],
        script=os.path.join(_SCRIPTS, "pon_mask_filter.py"),
    log:
        "logs/filtering/{sample}/pon_mask_filter.log",
    threads: max(2, workflow.cores // 4)
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=120,
    shell:
        """
        python {params.script} \
            --pon-fasta {params.pon_fasta} \
            --samtools-sif {params.samtools_sif} \
            --threads {threads} \
            {input.txt} > {output.txt} 2> {log}
        """
