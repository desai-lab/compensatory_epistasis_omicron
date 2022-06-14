#!/bin/bash
#SBATCH -p desai # Partition to submit to (comma separated)
#SBATCH -J epistasis # Job name
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-0:00 # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 100000 # Memory in MB
#SBATCH -o /n/holyscratch01/desai_lab/amoulana/omicron/Titeseq/epistasis_inference_test/20211214_ep_%A_%a.out
#SBATCH -e /n/holyscratch01/desai_lab/amoulana/omicron/Titeseq/epistasis_inference_test/slurm_output/20211214_ep_%A_%a.err

CURRDIR=$PWD

module load python/3.7.7-fasrc01
python infer_effects_H1_cv.py ${SLURM_ARRAY_TASK_ID}
python infer_effects_H1_final.py