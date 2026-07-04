#!/bin/bash
#SBATCH -A u2020001
#SBATCH -n 4
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=elena.vanzalen@umu.se
#SBATCH --mail-type=FAIL
#SBATCH -J "RunR"
set -ex

apptainer exec -B /mnt:/mnt -B /usr/local/lib/R/library R-4.4.3.sif  # provide your own R container, or call Rscript directly Rscript $@
