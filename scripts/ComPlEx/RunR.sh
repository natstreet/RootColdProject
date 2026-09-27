#!/bin/bash
#SBATCH -A PROJECT
#SBATCH -n 4
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=FAIL
#SBATCH -J "RunR"
# Runs an R script (its path and arguments are passed on) inside an R container; submitted by submit_ComPlEx_cold.sh.
# Provide your own R 4.4 container, or replace the apptainer line with: Rscript "$@"
set -ex
ulimit -c 0   # no core files if R crashes

apptainer exec -B /mnt:/mnt -B /usr/local/lib/R/library R-4.4.3.sif Rscript "$@"
