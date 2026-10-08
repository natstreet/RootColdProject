#!/bin/bash
# run_motif_tests.sh — the promoter motif test (D062 item 4, D063, D064) on the final core, for both windows and both
# modes (primary: genes without an annotated start codon excluded; sensitivity: included). Run in a working tree with
# data/motifs/ (motif_pairs.tsv.gz from motif_pairs.R, micro_out/ from motif_job_micro.sh and check_start_codons.py,
# jaspar2024_plants_clusters.tsv from jaspar_clusters.py) and results/core_stability_label.tsv (the 133 core, D066).
# Outputs: results/motifs/motif_{test,og_differences,pair_jaccard}_<win>_<mode>.*
set -euo pipefail
R=${RSCRIPT:-Rscript}; mkdir -p logs
for W in 1000 2000; do for M in primary sensitivity; do
  $R scripts/motifs/motif_conservation_test.R data/motifs/motif_pairs.tsv.gz data/motifs/micro_out data/motifs/jaspar2024_plants_clusters.tsv \
     results/core_stability_label.tsv $W $M results/motifs > logs/motif_${W}_${M}.log 2>&1
  echo "motif test $W $M: $(grep -c . results/motifs/motif_test_${W}_${M}.tsv) rows"
done; done
