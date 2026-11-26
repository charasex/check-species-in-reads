# check-species-in-reads

A Snakemake workflow to verify species identity and screen for cross-species contamination in paired-end sequencing reads. 

The pipeline supports two selectable approaches configured via `config/config.yaml`:
1. **Taxonomic classification (`kraken`)**: Builds a custom Kraken2 database from a reference FASTA (e.g., BOLD) and classifies paired-end reads.
2. **Targeted k-mer assembly & alignment (`kmer`)**: Extracts k-mers with `kmtricks`, retrieves matching reads via `back_to_sequences`, assembles them with `SPAdes`, and aligns contigs against a custom BOLD `BLASTn` database.

---

## Repository Structure

```text
check-species-in-reads/
├── snakefile                  # Main workflow entry point
├── config/
│   ├── config.yaml            # Workflow parameters and paths
│   └── samples.tsv            # Sample metadata sheet (species, sex, R1, R2)
└── workflow/
    ├── envs/
    │   ├── kmer.yaml          # Conda environment for the k-mer workflow
    │   └── kraken.yaml        # Conda environment for the Kraken2 workflow
    └── rules/
        ├── kmers.smk          # Rules for kmtricks, SPAdes, and BLASTn
        └── kraken.smk         # Rules for Kraken2 DB creation and classification
```

---

## Prerequisites

* [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
* [Snakemake](https://snakemake.readthedocs.io/) (tested with Snakemake v9.27.1)

All bioinformatics dependencies (`kraken2`, `kmtricks`, `back_to_sequences`, `spades`, `blast`) are automatically installed and managed per rule via `--use-conda`.

---

## Configuration

### 1. Sample Sheet (`samples.tsv`)
Create a tab-separated file (`samples.tsv`) at the folder `config/` of the repository linking each `(species_id, sex)` pair to its raw FASTQ files:

```tsv
species_id	sex	R1	R2
Leporinus_acutidens	M	Leporinus_acutidens_male-LACU-M_R1_001.fastq.gz	Leporinus_acutidens_male-LACU-M_R2_001.fastq.gz
Leporinus_acutidens	F	Leporinus_acutidens_female-LACU-F_R1_001.fastq.gz	Leporinus_acutidens_female-LACU-F_R2_001.fastq.gz
Leporinus_friderici	M	Leporinus_friderici_male-LFRI_M_R1_001.fastq.gz	Leporinus_friderici_male-LFRI_M_R2_001.fastq.gz
Leporinus_friderici	F	Leporinus_friderici_female-LFRI_F_R1_001.fastq.gz	Leporinus_friderici_female-LFRI_F_R2_001.fastq.gz
```

> **Note:** Input reads are expected to be located at `<wdir>/<species_id>/PoolSeq/MGX/<R1|R2>`.

### 2. Workflow Configuration (`config/config.yaml`)

```yaml
# Base directory containing species folders
wdir: "/path/to/data/raw"

# Workflow mode: "kraken" or "kmer"
tools: "kmer"

# Reference FASTA file (e.g., BOLD database)
BOLD: "/path/to/BOLD_reference.fasta"

# K-mer mode specific settings
kmerDir: "PoolSeq/MGX/check_contamination/kmers"
kmer_size: 21
```

---

## Usage

### 1. Dry-run (check DAG and commands without executing)
```bash
snakemake -n -p
```

### 2. Run the workflow
```bash
snakemake --use-conda --conda-frontend mamba --cores 10 --printshellcmds --rerun-incomplete
```

### 3. Switch between tools via the command line
You can override the `tools` parameter without editing `config/config.yaml`:

```bash
# Run the Kraken2 workflow
snakemake --use-conda --cores 10 --config tools=kraken

# Run the K-mer + SPAdes + BLASTn workflow
snakemake --use-conda --cores 10 --config tools=kmer
```

---

## Outputs

* **Databases** (built once in the repository root and reused across runs):
  * `CustomFishDB/` (Kraken2 index)
  * `CustomBlastBoldDB/` (BLASTn nucleotide index)
* **Per-sample results** (written inside `<wdir>/<species_id>/`):
  * **Kraken2 mode**: `PoolSeq/MGX/results/check_contamination/<species>_<sex>_kraken.report`
  * **K-mer mode**: `<kmerDir>/<species>_<sex>_blast_output.tsv` (tabular BLASTn `outfmt 6` results of assembled contigs against BOLD)

# tested tools
<img width="4200" height="2400" alt="detection_contamination_blast" src="https://github.com/user-attachments/assets/6ddfa009-9fbf-4f88-b3c1-529177ffa735" />
<img width="3600" height="2400" alt="kraken_taux_et_proportion" src="https://github.com/user-attachments/assets/912477bd-1de4-444e-b300-bd7d99557f63" />

