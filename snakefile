import pandas as pd

configfile: "config/config.yaml"

samples_df = pd.read_csv(config["sample"], sep="\t").set_index(["species_id", "sex"], drop=False)
WDIR = config["wdir"]

def get_fastqs(wildcards):
    row = samples_df.loc[(wildcards.species, wildcards.sex)]
    return {
        "r1": f"{WDIR}/{wildcards.species}/PoolSeq/MGX/{row['R1']}",
        "r2": f"{WDIR}/{wildcards.species}/PoolSeq/MGX/{row['R2']}",
    }

tool = config["tools"]

if tool == "kraken":
    include: "workflow/rules/kraken.smk"
    TARGETS = expand(
        WDIR + "/{species}/PoolSeq/MGX/check_contamination/{species}_{sex}_kraken.report",
        zip,
        species=samples_df["species_id"],
        sex=samples_df["sex"],
    )
elif tool == "kmer":
    include: "workflow/rules/kmers.smk"
    TARGETS = expand(
        WDIR + "/{species}/PoolSeq/MGX/check_contamination/kmers/{species}_{sex}_blast_output.tsv",
        zip,
        species=samples_df["species_id"],
        sex=samples_df["sex"],
    )
else:
    raise ValueError(f"unknow tools in config.yaml : '{tool}'. choose 'kraken' or 'kmer'.")



rule all:
    input:
        TARGETS