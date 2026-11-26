#!/bin/bash
#
# SLURM launcher. For PBS Pro or SGE, submit an equivalent wrapper with the
# scheduler's own directives and swap --profile for ./profiles/pbs or ./profiles/sge;
# see the "Running on another job scheduler" section of the README.
#
# Note the absence of --cores: it is set by the profile (cores: 64). Passing --cores N
# on the command line clamps EVERY rule's threads to N, even in cluster mode, so a
# --cores 1 here would silently undo the whole of set-threads.
#SBATCH --mail-user=ishaac.chninak@inrae.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name=kraken
#SBATCH --output=log/slurm-%A.out

# species=$1

# module load snakemake/8.9.0
# module load snakemake/8.27.1
# module load snakemake/9.4.0
# . /local/env/envconda.sh
# source activate snakemake
# module load conda
# source activate snakemake


# echo "Snakemake version"
# snakemake --version
# echo "SLURM version"
# srun --version
# echo "Conda version"
# conda --version

# git status

# mkdir --parents log

# echo "Running Snakemake pipeline for species $species..."

snakemake --configfile config/config.yaml --use-conda --conda-frontend mamba --printshellcmds --cores 10 --rerun-incomplete