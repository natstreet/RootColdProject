#!/bin/bash
# reproduce_figures_tables.sh — run the figure, table and validation producers, in order, in one working directory
# that holds (or links to) the analysis data: data/, DE/ (scripts/DE/differential_expression.R) and scripts/, plus the
# two STRING v12 A. thaliana files for Fig. S6b. Each producer's output stays in the working directory (figures and
# their value files), tables/ or results/; each producer's console output is kept in logs/.
#
# Usage: bash scripts/reproduce_figures_tables.sh <working dir> [Rscript for DESeq2 steps]
#   RSCRIPT: the Rscript to use (default: Rscript on the PATH). The second argument names an Rscript with DESeq2 for
#   tf_clique_directions.R (default: the same).
set -euo pipefail
W=$1; R=${RSCRIPT:-Rscript}; RD=${2:-$R}
cd "$W"; mkdir -p logs tables results
run() { local name=$1; shift; echo "== $name"; if "$@" > "logs/$name.log" 2>&1; then tail -3 "logs/$name.log"; else echo "FAILED: $name (logs/$name.log)"; tail -5 "logs/$name.log"; exit 1; fi; }
run genotype_colours      $R scripts/figures/genotype_colours.R
run TableS1               $R scripts/GO/make_TableS1.R
run TableS3               $R scripts/GO/make_TableS3.R
run TableS4               $R scripts/GO/make_TableS4.R
run Fig5_caption_cand     $R scripts/figures/fig5_caption_candidates.R
run Fig1                  $R scripts/figures/make_Fig1.R
run Fig2                  $R scripts/figures/make_Fig2.R
run Fig4                  $R scripts/figures/make_Fig4.R
run Fig5                  $R scripts/figures/make_Fig5.R Fig5_caption_candidates.tsv   # proposed terms; the final caption table is NRS's choice
run Fig6                  $R scripts/figures/make_Fig6.R
run FigS4                 $R scripts/figures/make_FigS4.R
run FigS5                 $R scripts/figures/make_FigS5.R
run FigS7                 $R scripts/figures/make_FigS7.R
run tf_clique_directions  $RD scripts/analysis/tf_clique_directions.R
run core_response         $R scripts/analysis/core_response_patterns.R
run core_modules          $R scripts/analysis/core_module_directions.R 3
run FigS2                 $R scripts/figures/make_FigS2.R OG0000166
run worked_example        $R scripts/analysis/worked_example_table.R OG0000166 results/worked_example_OG0000166.tsv
run cluster_profiles      $R scripts/analysis/cluster_profiles.R data/open results/cluster_profiles
run label_permutation     $R scripts/validation/conservation_label_permutation.R
run sizeaware_null        $R scripts/validation/conservation_sizeaware_null.R
run core_clusters         $R scripts/validation/make_core_clusters.R
run string_validation     $R scripts/validation/string_external_validation.R
run string_by_cluster     $R scripts/validation/string_by_cluster.R
run FigS6                 $R scripts/figures/make_FigS6.R
echo "all producers finished"
