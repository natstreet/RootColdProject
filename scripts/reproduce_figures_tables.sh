#!/bin/bash
# reproduce_figures_tables.sh — run the figure, table and validation producers, in order, in one working directory
# that holds (or links to) the analysis data: data/, DE/ (scripts/DE/differential_expression.R) and scripts/, plus the
# two STRING v12 A. thaliana files for Fig. S4b. Each producer's output stays in the working directory (figures and
# their value files), tables/ or results/; each producer's console output is kept in logs/.
#
# The orthology classification (D043) needs the N1 hierarchical orthogroups from the multiple-alignment gene trees and
# OrthoFinder's raw orthogroups in data/orthofinder/ (N1_msa.tsv, Orthogroups.tsv, msa_trees/).
#
# Sample stability (D052/D059 test as a label, D066) and its interpretation need: data/loo_run (the leave-one-sample-out
# ComPlEx runs: tasks.tsv, out/), data/stability_design/{matched_n,design} (outputs of stability_matched_n.R and
# design_stability.R, run on micro), and data/external_gene_sets (the published tables of external_gene_sets.R).
# The Fig. S1 worked example is drawn when WORKED_EXAMPLE names its orthogroup (OG0005042, D067). The companion-conifer analysis (A2,
# Fig. S5) refits DESeq2 on the dds objects of van Zalen et al. (2026, https://doi.org/10.64898/2026.08.03.742458) when COMPANION_DDS
# names their folder; otherwise it takes results made on a
# cluster from data/companion_a2/.
#
# Usage: bash scripts/reproduce_figures_tables.sh <working dir> [Rscript for DESeq2 steps]
#   RSCRIPT: the Rscript to use (default: Rscript on the PATH). The second argument names an Rscript with DESeq2 for
#   tf_clique_directions.R (default: the same). FIG4_CAPTIONS chooses the Fig. 4 captions: family_level (default;
#   fig4_captions_family_level.R on Table S3, D090), approved (Fig4_terms_approved.tsv, D008) or none.
set -euo pipefail
W=$1; R=${RSCRIPT:-Rscript}; RD=${2:-$R}
cd "$W"; mkdir -p logs tables results
run() { local name=$1; shift; echo "== $name"; if "$@" > "logs/$name.log" 2>&1; then tail -3 "logs/$name.log"; else echo "FAILED: $name (logs/$name.log)"; tail -5 "logs/$name.log"; exit 1; fi; }
run genotype_colours      $R scripts/figures/genotype_colours.R
run med_expression        $R scripts/expression/med_expression.R                         # R068 F02: per-timepoint VST medians (read by Figs 1, 4, S1 and the clusters)
run cluster_numbering     $R scripts/ComPlEx/manuscript_cluster_map.R --numbering            # D095: one numbering file (scripts/ComPlEx/cluster_numbering.tsv)
run cluster_origin        $R scripts/analysis/cluster_origin.R
run og_summary            $R scripts/DE/og_summary.R                                        # D101: deposited gene-to-orthogroup map (og_map.R)
run common_orthogroups    $R scripts/validation/common_orthogroups.R
run TableS1               $R scripts/GO/make_TableS1.R
run TableS2               $R scripts/GO/make_TableS2.R
run size_rule             $R scripts/ComPlEx/size_rule.R                                   # D026 even-odds rule table (D080); read by Table S3
run conifer_only          $R scripts/GO/conifer_only_go.R results/conifer_only            # R064 (D089): conifer-only co-expressologs, no DEG restriction
if [ -n "${HR_DOMAINS_INPUTS:-}" ]; then                                                  # R066: hmmscan needs Pfam-A and the proteomes ("<Pfam-A.hmm> <spruce aa> <pine aa>")
  run hr_domains          $R scripts/analysis/hr_family_domains.R results/hr_family_domains results/conifer_only $HR_DOMAINS_INPUTS
else echo "== hr_domains: HR_DOMAINS_INPUTS not set; results/hr_family_domains copied from data/hr_family_domains"; mkdir -p results/hr_family_domains; cp data/hr_family_domains/*.tsv results/hr_family_domains/; fi
run conifer_defence       $R scripts/deposit/conifer_defence_tables.R                      # R068: deposit tables for the D094 sentence
run TableS3_genotypes     $R scripts/GO/make_TableS3.R                                     # per genotype, gene level (supplementary comparison)
run TableS3_family_gt     $R scripts/GO/tableS3_family_level.R tables/TableS3_cliqueClusterGO_classic.tsv tables/family_level ${S3_NULL_DRAWS:-5000}
run TableS3               $R scripts/GO/tableS3_arabidopsis.R tables/arabidopsis            # D090: one test per cluster on the Col-0 annotation
run TableS3_install       bash -c 'mkdir -p tables/per_genotype && mv tables/TableS3_cliqueClusterGO_classic.tsv tables/per_genotype/TableS3_gene_level.tsv && mv tables/Table_S3.xlsx tables/per_genotype/Table_S3_gene_level.xlsx && cp tables/family_level/TableS3_family_level.tsv tables/per_genotype/TableS3_family_level.tsv && cp tables/family_level/Table_S3_family_level.xlsx tables/per_genotype/Table_S3_family_level.xlsx && cp tables/arabidopsis/TableS3_family_level.tsv tables/TableS3_cliqueClusterGO_classic.tsv && cp tables/arabidopsis/TableS3_gene_level.tsv tables/TableS3_gene_level.tsv && cp tables/arabidopsis/Table_S3.xlsx tables/Table_S3.xlsx && cp tables/arabidopsis/annotation_coverage.tsv tables/TableS3_annotation_coverage.tsv && cp tables/arabidopsis/annotation_coverage_summary.tsv tables/TableS3_annotation_coverage_summary.tsv && echo "Table S3 = Col-0 annotation, family level (D090); gene level tables/TableS3_gene_level.tsv; per-genotype runs in tables/per_genotype/"'
run Fig1                  $R scripts/figures/make_Fig1.R
run Fig2                  $R scripts/figures/make_Fig2.R
run Fig3                  $R scripts/figures/make_Fig3.R
run paralogue_agreement   $R scripts/analysis/paralogue_agreement.R results              # D083, drawn on Fig. 4
case "${FIG4_CAPTIONS:-family_level}" in                                                  # D090 captions (provisional until approved)
  none)     run Fig4      $R scripts/figures/make_Fig4.R none ;;
  approved) run Fig4_captions $R scripts/figures/fig4_captions_from_terms.R                # the approved terms (D008), matched to Table S3
            run Fig4      $R scripts/figures/make_Fig4.R Fig4_captions_approved.tsv ;;
  *)        run Fig4_captions $R scripts/figures/fig4_captions_family_level.R tables/TableS3_cliqueClusterGO_classic.tsv Fig4_captions_family_level.tsv 2
            run Fig4      $R scripts/figures/make_Fig4.R Fig4_captions_family_level.tsv ;;
esac
run size_flag             $R scripts/validation/orthogroup_size_flag.R
run hog_N1msa             $R scripts/analysis/hog_clique_consistency.R data/orthofinder/N1_msa.tsv results/hog_N1msa data/orthofinder/Orthogroups.tsv
run orthology_class       $R scripts/analysis/orthology_class.R results/hog_N1msa data/orthofinder/msa_trees
run loo_score             $R scripts/validation/loo_stability_score.R data/loo_run results/loo_stability   # D052/D059 test (gate first)
run stability_label       $R scripts/analysis/core_stability_label.R                                     # D066 label
run Fig6                  $R scripts/figures/make_Fig6.R                                   # TF cliques
run FigS2                 $R scripts/figures/make_FigS2.R
run FigS3                 $R scripts/figures/make_FigS3.R
run FigS6                 $R scripts/figures/make_FigS6.R                                  # CBF co-expressologs (R045 renumbering)
run tf_clique_directions  $RD scripts/analysis/tf_clique_directions.R
run core_response         $R scripts/analysis/core_response_patterns.R
run a1_rewiring           $R scripts/analysis/a1_rewiring_duplication.R                   # A1 (N159)
run aspen_degs            $R scripts/analysis/aspen_timepoint_degs.R                      # aspen DEGs per time point (R045)
if [ -n "${COMPANION_DDS:-}" ]; then                                                       # A2: DESeq2 refits of the companion conifer data (heavy; HPC)
  run a2_companion        $R scripts/companion/a2_tissue_stress_specificity.R "$COMPANION_DDS"
else echo "== a2_companion: COMPANION_DDS not set; results/a2_*.tsv copied from data/companion_a2 (made on micro)"; cp data/companion_a2/a2_*.tsv results/; fi
run FigS5                 $R scripts/figures/make_FigS5.R                                  # companion conifers (R045)
run orthology_tables      $R scripts/analysis/orthology_class_tables.R
run core_modules          $R scripts/analysis/core_module_directions.R 3
if [ -n "${WORKED_EXAMPLE:-}" ]; then
  run FigS1               $R scripts/figures/make_FigS1.R "$WORKED_EXAMPLE"
  run worked_example      $R scripts/analysis/worked_example_table.R "$WORKED_EXAMPLE" "results/worked_example_$WORKED_EXAMPLE.tsv"
else echo "== FigS1, worked_example: skipped (WORKED_EXAMPLE not set; choice pending)"; fi
run cluster_profiles      $R scripts/analysis/cluster_profiles.R data/open results/cluster_profiles
run Fig5                  $R scripts/figures/make_Fig5.R                                   # growth components (R043)
run label_permutation     $R scripts/validation/conservation_label_permutation.R
run sizeaware_null        $R scripts/validation/conservation_sizeaware_null.R
run core_clusters         $R scripts/validation/make_core_clusters.R
run string_validation     $R scripts/validation/string_external_validation.R
run string_by_cluster     $R scripts/validation/string_by_cluster.R
run FigS4                 $R scripts/figures/make_FigS4.R                                  # panel a: data/pseudo_family_null (R071, run on micro)
run loo_diagnostics       $R scripts/validation/loo_diagnostics.R data/loo_run results/loo_stability results/loo_stability/diagnostics
run kept_lost             $R scripts/validation/stability_kept_lost.R results/stability
run why_spruce            $R scripts/validation/stability_why_spruce.R results/stability
run timepoints            $R scripts/validation/stability_timepoint_attribution.R data/loo_run results/stability
run matched_n_summary     $R scripts/validation/stability_matched_n_summary.R data/stability_design/matched_n results/stability/matched_n_summary.tsv
run design_summary        bash -c "mkdir -p results/stability/design && cp data/stability_design/design/*.tsv results/stability/design/ && $R scripts/validation/design_stability_summary.R results/stability/design"
run string_stable         $R scripts/validation/string_external_validation.R results/stability/stable_core_clusters.tsv results/stability/stable
run string_sensitive      $R scripts/validation/string_external_validation.R results/stability/sensitive_core_clusters.tsv results/stability/sensitive
run external_sets         $R scripts/validation/external_gene_sets.R data/external_gene_sets data/external_gene_sets results/external_validation/sets
run external_test         $R scripts/validation/external_validation_test.R results/external_validation/sets/external_gene_sets.tsv results/stability/kept_lost_orthogroups.tsv results/external_validation
run menges2005            $R scripts/validation/menges2005_core_regulators.R data/external_gene_sets/tpj_2319_sm_tables1.xls results/stability/kept_lost_orthogroups.tsv results/external_validation
run worked_candidates     $R scripts/analysis/worked_example_candidates.R results/stability
run paralogue_switch      $R scripts/analysis/paralogue_switching.R                       # D062 item 2
run paralogue_artefact    $R scripts/analysis/paralogue_switching_artefact.R              # needs data/paralogue_switching (Araport11 extracts)
run motif_tests           bash scripts/motifs/run_motif_tests.sh                           # D062 item 4 (data/motifs)
run cold_motif_check     $R scripts/motifs/cold_motif_check.R data/cbf_check/promoters scripts/motifs/cold_panel_motifs.tsv results/external_validation/sets/external_gene_sets.tsv data/cbf_check/orthogroup_core_status.tsv results/cbf_check
run ga_biosynthesis       $R scripts/analysis/ga_biosynthesis_check.R                     # D049 GA paragraph
run manuscript_values     $R scripts/analysis/manuscript_values.R
echo "all producers finished"
