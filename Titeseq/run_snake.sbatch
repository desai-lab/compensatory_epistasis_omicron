#!/bin/bash
#
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J RBD_abecgivdhjx # Job name
#SBATCH --ntasks 48
#SBATCH -t 0-12:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 64000 # Memory in MB
#SBATCH -o /n/desai_lab/users/tdupic/RBDabecgivdjhx_snakemake_%A_%a.out
#SBATCH -e /n/desai_lab/users/tdupic/RBDabecgivdjhx_snakemake_%A_%a.err

module load python/3.7.7-fasrc01
source activate omicron

snakemake -j 48
