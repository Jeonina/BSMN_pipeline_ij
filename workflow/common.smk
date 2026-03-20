import pandas as pd

# --- Load config and samples ------------------------------------------------

configfile: "config/config.yaml"

samples_df = pd.read_csv(config["samples"], sep="\t", dtype=str)
samples_df = samples_df.set_index(["sample_id", "readgroup"], drop=False)

SAMPLES = samples_df["sample_id"].unique().tolist()

REF = config["ref"]["fasta"]

# --- Helper functions --------------------------------------------------------

def get_fastqs(wildcards):
    """Return fq1, fq2 for a given sample + readgroup."""
    row = samples_df.loc[(wildcards.sample, wildcards.rg)]
    return {"fq1": row["fq1"], "fq2": row["fq2"]}


def get_readgroups(sample):
    """Return list of readgroups for a sample."""
    return samples_df.loc[sample, "readgroup"].unique().tolist()


def get_sorted_bams(wildcards):
    """Return all sorted BAMs for a sample (one per readgroup)."""
    rgs = get_readgroups(wildcards.sample)
    return expand(
        "results/mapping/{sample}/align/{sample}.{rg}.sorted.bam",
        sample=wildcards.sample, rg=rgs,
    )
