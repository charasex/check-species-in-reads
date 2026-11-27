import pandas as pd

configfile: "config/config.yaml"

samples_df = pd.read_csv(config["sample"], sep="\t").set_index(["species_id", "sex"], drop=False)
# path to the read files and output directory
WDIR = config["wdir"]
OUTDIR = config["outdir"]

# function to get the fastq files for a given species and sex
def get_fastqs(wildcards):
    row = samples_df.loc[(wildcards.species, wildcards.sex)]
    return {
        "r1": f"{WDIR}/{wildcards.species}/PoolSeq/MGX/{row['R1']}",
        "r2": f"{WDIR}/{wildcards.species}/PoolSeq/MGX/{row['R2']}",
    }

# Checksums

CHECKSUM_DIR = OUTDIR + "/{species}/checksums/{sex}"

INITIAL_CHECKSUMS = CHECKSUM_DIR + "/initial.sha256"
VERIFIED_CHECKSUMS = CHECKSUM_DIR + "/integrity_verified.txt"


tool = config["tools"]

# choose the workflow to run based on the tool specified in the config file
if tool == "both":
    include: "workflow/rules/kraken.smk"
    include: "workflow/rules/kmers.smk"
    include: "workflow/rules/plot_result.smk"
    TARGETS = (
        expand(OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_kraken.report", 
            zip, 
            species=samples_df["species_id"], 
            sex=samples_df["sex"],
        )
        + expand(OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_blast_output.tsv", 
            zip,
            species=samples_df["species_id"], 
            sex=samples_df["sex"],
        )
        + expand(OUTDIR + "/{species}/PoolSeq_check/plots/{species}_{sex}_blast.png",
            zip,
            species=samples_df["species_id"], 
            sex=samples_df["sex"],
        )
        + expand(OUTDIR + "/{species}/PoolSeq_check/plots/{species}_{sex}_kraken.png",
            zip,
            species=samples_df["species_id"], 
            sex=samples_df["sex"],
        )
    )

elif tool == "kraken":
    include: "workflow/rules/kraken.smk"
    include: "workflow/rules/plot_result.smk"
    TARGETS = (
        expand(OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_kraken.report",
            zip,
            species=samples_df["species_id"],
            sex=samples_df["sex"],
        )
        + expand(OUTDIR + "/{species}/PoolSeq_check/plots/{species}_{sex}_kraken.png",
            zip,
            species=samples_df["species_id"], 
            sex=samples_df["sex"],
        )
    )
elif tool == "kmer":
    include: "workflow/rules/kmers.smk"
    include: "workflow/rules/plot_result.smk"
    TARGETS = (
        expand(OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_blast_output.tsv",
            zip,
            species=samples_df["species_id"],
            sex=samples_df["sex"],
        )
        + expand(OUTDIR + "/{species}/PoolSeq_check/plots/{species}_{sex}_blast.png",
            zip,
            species=samples_df["species_id"], 
            sex=samples_df["sex"],
        )
    )
else:
    raise ValueError(f"unknow tools in config.yaml : '{tool}'. choose 'kraken' or 'kmer' or 'both'")


include: "workflow/rules/checksum.smk"

rule all:
    input:
        TARGETS,
        expand(
            VERIFIED_CHECKSUMS,
            zip,
            species=samples_df["species_id"],
            sex=samples_df["sex"],
        )