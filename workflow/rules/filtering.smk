# =============================================================================
# Filtering rules: Mutect2-filtered VCF → final somatic SNV list
#
# Default cascade — six sequential filters per PIPELINE.md:
#   1. accessibility_filter   — 1KG strict mask (keeps only mask base == 'P')
#   2. germline_filter        — removes gnomAD AF > 0.001 known germline variants
#   3. vaf_filter             — binomial test p < 1e-6 AND alt_count >= 5
#   4. cnvnator_filter        — BSMN D-step: drop CNV-region calls (CN >= 2.5)
#   5. mosaicforecast_filter  — BSMN E-step: MosaicForecast RF mosaic prediction
#   6. pon_mask_filter        — IUPAC FASTA-based panel-of-normals mask
#
# BSMN E-step is mayo OR MosaicForecast — two alternatives, NOT a chain. The
# original BSMN runs each as its own branch (jobs/submit_filtering_jobs.py).
# We default to MosaicForecast only. The mayo_filter rule + scripts + tests are
# retained but NOT wired into the default cascade; to enable mayo instead, point
# mosaicforecast_filter.input.txt back at {sample}.mayo_filtered.txt.
#
# cnvnator_root (prep) builds the per-sample read-depth ROOT consumed by
# cnvnator_filter. Resources are dynamically allocated (config/resolved_params.yaml).
#
# DAG (per sample):
#   results/calling/{sample}/{sample}.filtered.vcf.gz
#       → accessibility_filter   → {sample}.accessible.txt (temp)
#       → germline_filter        → {sample}.germline_filtered.txt (temp)
#       → vaf_filter             → {sample}.vaf_filtered.txt (temp)
#       → cnvnator_filter        → {sample}.cnvnator_filtered.txt (temp)
#           (uses cnvnator/{sample}.root from rule cnvnator_root)
#       → mosaicforecast_filter  → {sample}.mosaicforecast_filtered.txt (temp)
#       → pon_mask_filter        → {sample}.final.txt  ← final output
#   (mayo_filter is defined but unwired by default — see note above)
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


rule cnvnator_root:
    """
    Build the per-sample CNVnator read-depth ROOT (BSMN A.CNVnator_mk_root.sh).

    CNVnator reads alignments via htslib; CRAM decoding is done reliably by
    samtools (-T ref) into a temp BAM restricted to the analysis chromosomes,
    which is then fed to `-tree`. v0.4.1 builds the GC histogram from a single
    reference via `-his ... -fasta` (no per-chromosome split needed).
    """
    input:
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        root="results/filtering/{sample}/cnvnator/{sample}.root",
    params:
        ref=REF,
        chrom=" ".join(CHROMOSOMES),
        binsize=_filtering.get("cnvnator", {}).get("binsize", 100),
        outdir="results/filtering/{sample}/cnvnator",
        tmpbam="results/filtering/{sample}/cnvnator/{sample}.tmp.bam",
        cnvnator_sif=CONTAINERS["cnvnator"]["sif"],
        samtools_sif=CONTAINERS["samtools"]["sif"],
    log:
        "logs/filtering/{sample}/cnvnator_root.log",
    threads: max(4, workflow.cores // 2)
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024 * 2,
        runtime=1440,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[cnvnator_root] START $(date -Iseconds)  sample={wildcards.sample}"
        echo "[cnvnator_root] chrom={params.chrom}  binsize={params.binsize}"
        echo "================================================================"
        mkdir -p {params.outdir}
        rm -f {output.root}

        # Reliable CRAM decode (samtools -T ref) → temp BAM on analysis chromosomes
        apptainer exec {params.samtools_sif} \
            samtools view -b -@ {threads} -T {params.ref} {input.cram} {params.chrom} \
            > {params.tmpbam}
        apptainer exec {params.samtools_sif} samtools index {params.tmpbam}

        # CNVnator ROOT: tree → his(-fasta) → stat → partition → call
        apptainer exec {params.cnvnator_sif} cnvnator -root {output.root} -chrom {params.chrom} -tree {params.tmpbam} -lite
        apptainer exec {params.cnvnator_sif} cnvnator -root {output.root} -chrom {params.chrom} -his {params.binsize} -fasta {params.ref}
        apptainer exec {params.cnvnator_sif} cnvnator -root {output.root} -chrom {params.chrom} -stat {params.binsize}
        apptainer exec {params.cnvnator_sif} cnvnator -root {output.root} -chrom {params.chrom} -partition {params.binsize}
        apptainer exec {params.cnvnator_sif} cnvnator -root {output.root} -chrom {params.chrom} -call {params.binsize} > {params.outdir}/{wildcards.sample}.cnvcall

        rm -f {params.tmpbam} {params.tmpbam}.bai
        echo "[cnvnator_root] END $(date -Iseconds)"
        echo "================================================================"
        """


rule cnvnator_filter:
    """
    CNVnator genotype filter — BSMN D.CNVnator_genotype_filter.sh step.
    Genotypes a +/-1kb window around each VAF-passing candidate and drops those
    whose estimated copy number is >= cn_threshold (2.5) — i.e. in a duplicated
    region where low apparent VAF is a copy-number / paralog artifact.
    """
    input:
        txt="results/filtering/{sample}/{sample}.vaf_filtered.txt",
        root="results/filtering/{sample}/cnvnator/{sample}.root",
    output:
        txt=temp("results/filtering/{sample}/{sample}.cnvnator_filtered.txt"),
    params:
        binsize=_filtering.get("cnvnator", {}).get("binsize", 100),
        cn_threshold=_filtering.get("cnvnator", {}).get("cn_threshold", 2.5),
        window=_filtering.get("cnvnator", {}).get("window", 1000),
        cnvnator_sif=CONTAINERS["cnvnator"]["sif"],
        script=os.path.join(_SCRIPTS, "cnvnator_filter.py"),
    log:
        "logs/filtering/{sample}/cnvnator_filter.log",
    threads: 1
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024,
        runtime=240,
    shell:
        """
        exec >> {log} 2>&1
        echo "[cnvnator_filter] START $(date -Iseconds)  sample={wildcards.sample}"
        python {params.script} \
            --candidates {input.txt} \
            --root {input.root} \
            --cnvnator-sif {params.cnvnator_sif} \
            --binsize {params.binsize} \
            --cn-threshold {params.cn_threshold} \
            --window {params.window} \
            > {output.txt}
        echo "[cnvnator_filter] kept=$(wc -l < {output.txt})"
        echo "[cnvnator_filter] END $(date -Iseconds)"
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
        txt="results/filtering/{sample}/{sample}.cnvnator_filtered.txt",
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


rule mosaicforecast_filter:
    """
    MosaicForecast filter — BSMN E.MosaicForecast.sh step.

    Runs MosaicForecast read-level feature extraction + trained-RF prediction
    (via the yanmei/mosaicforecast Apptainer image) on the CNVnator-passing
    candidates and keeps only those predicted 'mosaic'. The MF scripts and the
    k24 mappability bigwig live inside the image; the RF model is on the host
    (cloned MosaicForecast repo). Apptainer auto-binds the project working dir,
    so the reference/CRAM/model/temp paths resolve unchanged inside the
    container. min_prob=0 keeps all mosaic calls (BSMN OUT_ALL); set 0.6 for
    the high-confidence set (OUT_HC).
    """
    input:
        txt="results/filtering/{sample}/{sample}.cnvnator_filtered.txt",
        cram="results/mapping/{sample}/{sample}.cram",
        crai="results/mapping/{sample}/{sample}.cram.crai",
    output:
        txt=temp("results/filtering/{sample}/{sample}.mosaicforecast_filtered.txt"),
    params:
        ref=REF,
        bam_dir="results/mapping/{sample}",
        workdir="results/filtering/{sample}/mf",
        model=_filtering.get("mosaicforecast", {}).get(
            "model", "resources/MosaicForecast/models_trained/250xRFmodel_addRMSK_Refine.rds"
        ),
        mode=_filtering.get("mosaicforecast", {}).get("mode", "Refine"),
        min_prob=_filtering.get("mosaicforecast", {}).get("min_prob", 0.0),
        timeout=_filtering.get("mosaicforecast", {}).get("timeout", 300),
        retries=_filtering.get("mosaicforecast", {}).get("retries", 5),
        mf_sif=CONTAINERS["mosaicforecast"]["sif"],
        script=os.path.join(_SCRIPTS, "mosaicforecast_filter.py"),
    log:
        "logs/filtering/{sample}/mosaicforecast_filter.log",
    threads: max(2, workflow.cores // 4)
    resources:
        mem_mb=lambda wildcards: _gatk_mem_gb * 1024 * 2,
        runtime=720,
    shell:
        """
        exec >> {log} 2>&1
        echo "================================================================"
        echo "[mosaicforecast_filter] START $(date -Iseconds)"
        echo "[mosaicforecast_filter] sample={wildcards.sample}  model={params.model}"
        echo "================================================================"

        python {params.script} \
            --candidates {input.txt} \
            --sample {wildcards.sample} \
            --bam-dir {params.bam_dir} \
            --fmt cram \
            --ref {params.ref} \
            --model {params.model} \
            --mf-sif {params.mf_sif} \
            --workdir {params.workdir} \
            --mode {params.mode} \
            --threads {threads} \
            --timeout {params.timeout} \
            --retries {params.retries} \
            --min-prob {params.min_prob} \
            > {output.txt}

        _kept=$(wc -l < {output.txt})
        echo "[mosaicforecast_filter] variants_kept=$_kept"
        echo "[mosaicforecast_filter] END $(date -Iseconds)"
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
        txt="results/filtering/{sample}/{sample}.mosaicforecast_filtered.txt",
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
