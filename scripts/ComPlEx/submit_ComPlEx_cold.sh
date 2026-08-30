#!/bin/bash
# Submit one ComPlEx pairwise co-expressolog comparison (set the species pair in ComPlEx.R)
# Run from: ./
set -eux

proj=PROJECT
mail=you@example.org
rscript=scripts/ComPlEx/ComPlEx.R
out=data/ComPlEx

if [ ! -d "$out" ]; then mkdir -p "$out"; fi

sbatch -A "$proj" -t 2-00:00:00 --mem=200G -n 4 \
  --mail-user="$mail" --mail-type=END,FAIL \
  -e "$out/ComPlEx.err" -o "$out/ComPlEx.out" \
  -J ComPlEx \
  scripts/ComPlEx/RunR.sh "$rscript"

echo "Submitted ComPlEx job. Output: $out"
