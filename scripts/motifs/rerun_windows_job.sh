#!/bin/bash
#SBATCH -J motif_win2
#SBATCH -p main
#SBATCH --exclude=R2D2,pando,kalkyl
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o logs/motif_win2_%j.out
# Re-runs promoter_windows.R for the five genotypes (adds the first3 column for check_start_codons.py, D064), checks that
# the window BEDs are byte-identical to the previous run, then writes the per-gene start-codon flags.
ulimit -c 0
set -euo pipefail
W=/mnt/ada/projects/spruce/nstreet/motifs_20261008; cd $W
export PIXI_HOME=$W/.pixi_home PIXI_CACHE_DIR=$W/.pixi_cache PATH=$W/.pixi/envs/default/bin:$PATH
mkdir -p out/prev_windows; [ -s out/prev_windows/col0_1000.bed ] || cp out/*_windows.tsv out/*_1000.bed out/*_2000.bed out/prev_windows/
while IFS=$'\t' read -r g gff fa pfx rep hard; do
  /usr/bin/time -f %e_s_%M_kB Rscript scripts/promoter_windows.R "$gff" inputs/genes_$g.txt "$fa.fai" "$pfx" "$rep" out/$g
  for w in 1000 2000; do cmp out/${g}_$w.bed out/prev_windows/${g}_$w.bed; done; echo "$g windows unchanged"
done < tmp/genomes.tsv
python3 scripts/check_start_codons.py > out/start_codons.tsv
awk -F'\t' 'NR>1{n[$1]++; if($6=="TRUE") a[$1]++} END{for(g in n) print g, n[g], a[g]+0}' out/start_codons.tsv
