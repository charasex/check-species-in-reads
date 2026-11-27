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
├── run_snakemake.sh           # Bash script to run in background the pipeline
├── config/
│   ├── config.yaml            # Workflow parameters and paths
│   └── samples.tsv            # Sample metadata sheet (species, sex, R1, R2)
└── workflow/
    ├── envs/
    │   ├── kmer.yaml          # Conda environment for the k-mer workflow
    │   ├── kraken.yaml        # Conda environment for the Kraken2 workflow
    │   └── r_plots.yaml       # Conda environment with package R for plot
    ├── rules/
    │    ├── checksum.smk       # Rules for check the integrity of input
    │    ├── kmers.smk          # Rules for kmtricks, SPAdes, and BLASTn
    │    ├── kraken.smk         # Rules for Kraken2 DB creation and classification
    │    └── plot_result.smk    # Rules for plotting result of output of kraken and blast
    └── scripts/
        ├── plot_result_blast.R    # Script for generate plot
        └── plot_result_kraken.R   # Script for generate plot
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
> 
> You need to modify the path of your input in config file and also adapt the path above

### 2. Workflow Configuration (`config/config.yaml`)

```yaml
# Base directory containing species folders
wdir: "/path/to/data/raw"
outdir: "/paht/to/result/output"

# Workflow mode: "kraken" or "kmer"
tools: "both" #or use specific kraken or kmer

# Reference FASTA file (e.g., BOLD database)
BOLD: "/path/to/BOLD_reference.fasta" # for Kraken
BOLD_blast: "/path/to/BOLD_reference_blast.fasta" # 
# K-mer mode specific settings
kmer_size: 21
```

---

## Usage

### 1. Dry-run (check DAG, commands and path without executing)
```bash
snakemake -n -p
```

### 2. Run the workflow
```bash
snakemake --use-conda --conda-frontend mamba --cores 10 --printshellcmds

# or use bash script for that
sbatch run_snakemake.sh
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
  * **Kraken2 mode**: `<OUTDIR>/PoolSeq_check/<species>_<sex>_kraken.report`
  * **K-mer mode**: `<OUTDIR>/PoolSeq_check/<species>_<sex>_blast_output.tsv` (tabular BLASTn `outfmt 6` results of assembled contigs against BOLD)
* **Per-sample Plot**
  * **For kraken2**: `<OUTDIR>/PoolSeq_check/plots/<species>_<sex>_kraken.png`
  * **For K-mer mode**: `<OUTDIR>/PoolSeq_check/plots/<species>_<sex>_blast.png`

## tested tools
We have tested the confidence of our methods to see if we can trust the result. In this case, we create 4 sample 0/30 (no contamination), 1/30 (contamination 1 by 30), 5/30 (contamination 5 by 30) and 15/30 (contamination 15 by 30).
<img width="4200" height="2400" alt="detection_contamination_blast" src="https://github.com/user-attachments/assets/6ddfa009-9fbf-4f88-b3c1-529177ffa735" />
<img width="3600" height="2400" alt="kraken_taux_et_proportion" src="https://github.com/user-attachments/assets/912477bd-1de4-444e-b300-bd7d99557f63" />

