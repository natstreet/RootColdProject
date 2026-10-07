# Reproducing the figures and tables

The producer scripts in `scripts/` regenerate the manuscript's figures and
supplementary tables from the processed data deposited on figshare
(doi:10.17044/scilifelab.32747586).

## Setup

1. Download `data.tar.gz` from the figshare deposit and unpack it in the repository
   root so that the data sits under `data/`:

   ```
   tar xzf data.tar.gz          # creates ./data/...
   ```

2. Run any producer from the repository root, e.g.:

   ```
   Rscript scripts/GO/make_TableS4.R
   ```

Outputs are written to the working directory. R package requirements: `DESeq2`,
`matrixStats`, `topGO`, `GO.db`, `ggplot2`, `ggtext`, `corrplot`, `dplyr`, `tidyr`,
`pheatmap`, `patchwork`, `cowplot`, `ggplotify`, `VennDiagram`, `UpSetR`, `readr`, `readxl`, `writexl`,
`data.table` and `curl`.

Some small derived inputs are provided in this repository rather than in the figshare deposit:
the median-expression tables (`expression/med_expression/`), the gene→orthogroup maps
(`ComPlEx/gene2OG/`), the per-species GO backgrounds (`annotation/bg_*.rda`), the Arabidopsis gene
aliases (`annotation/gene_aliases_20140331.txt`), the per-orthogroup DEG summary
(`DEGs/og_summary.RData`, from `scripts/DE/og_summary.R`), and the per-timepoint DEG lists (`DEGs/per_timepoint_DEG_lists/`). Unpacking the deposit alongside them gives a complete `data/` tree.

## Definitions used throughout

- **DEG:** adjusted P ≤ 0.05 and |log2 fold change| ≥ 1 (a two-fold change). The per-timepoint
  lists (`…_p0.05_lfc0.txt`, provided in this repository) and the deposited per-species union sets
  (`DEGs_*.RData`) are the DEGs used throughout the paper.
- **GO enrichment:** topGO, classic Fisher's exact test, biological process,
  Benjamini–Hochberg adjustment, padj < 0.05 (Tables S1, S3 and S4).

## Producer map

| Output | Script | Key inputs (under `data/`) |
|---|---|---|
| Genotype colours (every figure that colours by genotype; the deposit explorer) | `scripts/figures/genotype_colours.R` (writes `genotype_colours.json`); `check_genotype_colours.R` fails if a producer defines its own | none |
| Fig. 1 | `scripts/figures/make_Fig1.R` | `expression/med_expression/`, `ComPlEx/gene2OG/`, `annotation/bg_*`, `DEGs/`; panel values in `results/Figure1_panel_values.tsv` |
| Fig. 2 | `scripts/figures/make_Fig2.R` | `DEGs/per_timepoint_DEG_lists/`, `DEGs/log2FoldChange_*.RData` |
| Fig. 4 | `scripts/figures/make_Fig4.R` | `ComPlEx/co_expressologs.RData`; panel b from the Table S3 output (`tables/TableS3_cladeGO_classic_p10.tsv`); Venn regions in `results/Figure4_venn_regions.tsv`, panel b values in `results/Figure4b_values.tsv` |
| Fig. 5 | `scripts/figures/make_Fig5.R` | `ComPlEx/cliques/{conserved_genes_COMPLETE_filtered.RDS,clique_HMcluster_list.rda}`, `expression/med_expression/`; the Table S4 output and `Fig5_captions_approved.tsv` (block captions by GO ID), written by `scripts/figures/fig5_captions_from_terms.R` from the approved terms in `scripts/figures/Fig5_terms_approved.tsv` (D008) |
| Fig. 6 | `scripts/figures/make_Fig6.R` | `annotation/gene_aliases_20140331.txt`, `ComPlEx/cliques/{clique_genes_filterable_COMPLETE.RDS,conserved_genes_COMPLETE_filtered.RDS}`, `expression/`, `DEGs/`; `results/orthogroup_size_flag.tsv` and `results/orthology_class.tsv` (panel labels) |
| Fig. S2 | `scripts/figures/make_FigS2.R <orthogroup>` | `ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS`, `DE/de_results.tsv`, `expression/`, `ComPlEx/orthologs-*.RData`, `annotation/{Athal_go_ids.tsv,gene_aliases_20140331.txt}`; `tables/core_conserved_neighbours.tsv.gz` from `scripts/ComPlEx/core_neighbourhoods.R` (15 pairs) and `scripts/analysis/core_module_directions.R`; `results/orthology_class.tsv`. The worked example is NOP56, OG0004969 (D044) |
| Fig. S4 | `scripts/figures/make_FigS4.R` | `ComPlEx/co_expressologs.RData`, `DEGs/DEGs_*.RData` |
| Fig. S5 | `scripts/figures/make_FigS5.R` | `ComPlEx/co_expressologs.RData`; per-comparison counts in `results/Figure_S5_values.tsv` |
| Fig. S6 | `scripts/figures/make_FigS6.R` | outputs of `scripts/validation/conservation_label_permutation.R` (panel a) and `make_core_clusters.R` + `string_external_validation.R` (panel b; **STRING download, see note**) |
| Fig. S7 | `scripts/figures/make_FigS7.R` | `ComPlEx/co_expressologs.RData`, `expression/` |
| Table S1 (`Table_S1.xlsx`) | `scripts/GO/make_TableS1.R` (with `go_table_utils.R`) | `DEGs/og_summary.RData`, `DEGs/DEGs_col0.RData`, `*/Orthogroups_20240613.tsv`, `annotation/Athal_go_ids.tsv`, `annotation/bg_col0_ATC.rda` |
| Table S3 (`Table_S3.xlsx`) | `scripts/GO/make_TableS3.R` (with `go_table_utils.R`) | `ComPlEx/co_expressologs.RData`, `annotation/{go_ids,bg_*}`, `DEGs/DEGs_*.RData` |
| Table S4 (`Table_S4.xlsx`) | `scripts/GO/make_TableS4.R` (with `go_table_utils.R`) | `ComPlEx/cliques/{clique_HMcluster_list.rda,clique_genes_filterable_COMPLETE.RDS}`, `annotation/{go_ids,bg_*}` |
| DEGs (Fig. 1 / Fig. 2) | `scripts/DE/differential_expression.R` | `dds/dds_{col0_soil,ost0_soil,Pt_new,Bp_new,Pa,Ps}.rda` |
| GO backgrounds | `scripts/annotation/make_backgrounds.R` | `expression/<sp>_root_cold_stress_expression.txt` |
| Species networks | `scripts/networks/export_network_inputs.R`, then `scripts/networks/build_coexpression_network.R` per species (method: `docs/coexpression_networks.md`) | `dds/` |
| TF clique responses (Discussion) | `scripts/analysis/tf_clique_directions.R` | `dds/`, `expression/`, `ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS`, `annotation/gene_aliases_20140331.txt`; output `results/tf_clique_directions.tsv` |
| Orthogroups split into several HOGs (Methods) | `scripts/analysis/orthogroup_hog_split.R` | `ComPlEx/Orthogroups_20240613.tsv`, `ComPlEx/Phylogenetic_Hierarchical_Orthogroups_N0.tsv.gz` (in this repository; N0.tsv of the same OrthoFinder run, from github.com/tallgran/conifer-JasonHill); output `results/orthogroup_hog_split.tsv` |
| Complete cliques per core orthogroup (Methods) | `scripts/analysis/complete_cliques_per_orthogroup.R` | `ComPlEx/cliques/{clique_genes_filterable_COMPLETE.RDS,conserved_genes_COMPLETE_filtered.RDS}`; outputs `results/complete_cliques_per_orthogroup.tsv`, `results/complete_cliques_summary.tsv` |
| Cluster profiles (Results, clusters 1–2) | `scripts/analysis/cluster_profiles.R` | `open/{gene_index,conserved_core,orthogroup_annotation}.tsv` — **deposit v4 or later**; outputs in `results/cluster_profiles/` |
| Common DEG orthogroups: count per gene-to-orthogroup map, and share with ≥ 20 genes (Results) | `scripts/validation/common_orthogroups.R` | `DEGs/og_summary.RData`, `DEGs/DEGs_*.RData`, `*/Orthogroups_20240613.tsv`, `ComPlEx/gene2OG/`, `open/gene_orthogroup_map.tsv` (deposit v4 or later; optional); output `results/common_orthogroups.tsv` |
| Col-0 DEGs that are also Ost-0 DEGs (Results) | `scripts/validation/ecotype_deg_overlap.R` | `DEGs/per_timepoint_DEG_lists/{col0,ost0}/`, `DEGs/DEGs_{col0,ost0}.RData`; output `results/ecotype_deg_overlap.tsv` |
| Co-expressolog comparisons per orthogroup (Results; run with `0.1 OG0000037 OG0001369`) | `scripts/validation/orthogroup_coexpressolog_comparisons.R` | `ComPlEx/co_expressologs.RData`, `open/orthogroup_annotation.tsv` — **deposit v4 or later**; output `results/orthogroup_coexpressolog_comparisons.tsv` |
| Family-size-aware null for the conserved core (Results) | `scripts/validation/conservation_sizeaware_null.R` | `open/co_expressologs.tsv` — **deposit v4 or later**, `*/Orthogroups_20240613.tsv`; outputs `results/conservation_sizeaware_{null.csv,summary.tsv}` |
| ENA runs, sequencing depth and platform per sample (Methods; deposit `samples.tsv`) | `scripts/deposit/sample_ena_runs.R` | `open/samples.tsv`, `quantification/{sample_name_map,lane_map}.tsv` — **deposit v4 or later**; ENA portal API (network); outputs `results/{sample_ena_runs,ena_run_metadata,sequencing_depth}.tsv` |
| Per-gene DE results (fold change, padj against no change, DEG call) | `scripts/DE/differential_expression.R` | as for DEGs; output `DE/de_results.tsv` |
| Deposit GO tables (`open/go_enrichment_*.tsv`) | `scripts/deposit/export_go_tables.R` | the Table S1, S3 and S4 workbooks |
| Deposit manuscript cluster numbers (`open/conserved_core.tsv`, `open/gene_index.tsv`) | `scripts/deposit/manuscript_clusters.R` | `MAP` in `scripts/GO/make_TableS4.R` |
| Deposit DEG calls (`open/differential_expression.tsv` `padj_lfc0`, `is_deg`; `open/gene_index.tsv` cold-response columns) | `scripts/deposit/de_calls.R` | `DE/de_results.tsv`, `DEGs/DEGs_*.RData` |
| Deposit ENA runs per sample (`open/samples.tsv`) | `scripts/deposit/sample_runs_columns.R` | `results/sample_ena_runs.tsv` |
| Deposit open tables from the analysis objects (`open/{co_expressologs,conserved_core,pairwise_comparisons,orthogroup_annotation,gene_orthogroup_map,conserved_vs_specific_genes,differential_expression,transcription_factors}.tsv`) | `scripts/deposit/export_open_formats.R` (also shipped in the deposit) | `ComPlEx/`, `DEGs/`, `TranscriptionFactors/`, `orthologs/orthologs.parquet` |
| Deposit size rule and orthology class (`open/conserved_core.tsv` last columns; `open/conserved_core_flagged.tsv`) | `scripts/deposit/core_size_columns.R` | `results/orthogroup_size_rule.tsv`, `results/orthology_class.tsv`, `ComPlEx/cliques/conserved_genes_COMPLETE_filtered_flagged.RDS` |
| Deposit sample table, design columns (`open/samples.tsv`) | `scripts/deposit/samples_table.R` | the expression-matrix column names |
| Deposit core network for Cytoscape (`open/conserved_core_network.sif`) | `scripts/deposit/core_network_sif.R` | `open/co_expressologs.tsv`, `open/conserved_core.tsv` |
| Deposit gene index (`open/gene_index.tsv`) | `scripts/deposit/gene_index.R` | the other open tables |
| Deposit Fig. 2 source values (`figure_source_values/Fig2_perTimepoint_lfc0.csv`) | `scripts/figures/make_Fig2.R` (`figures/Fig2_perTimepoint_lfc0.csv`) | as Fig. 2 |

## Upstream: how the deposited analysis objects were made

The producers above start from the deposited data. These scripts made the deposited analysis objects from the
expression data, and are needed only to rebuild those objects:

| Deposited object | Script |
|---|---|
| `expression/<species>_root_cold_stress_expression.txt` (variance-stabilised matrices) | `scripts/ComPlEx/Network_Generation_setup.R` |
| `ComPlEx/PairwiseComparisons/` (the 15 pairwise co-expressolog comparisons) | `scripts/ComPlEx/ComPlEx.R`, one run per species pair, submitted with `scripts/ComPlEx/submit_ComPlEx_cold.sh` (which uses `RunR.sh`) |
| `ComPlEx/co_expressologs.RData`, `ComPlEx/all_expressed_genes.RData` | `scripts/ComPlEx/Genelists.Rmd` |
| `ComPlEx/ones_and_zeros_all.RData` (per orthogroup and comparison, whether ComPlEx tested a gene pair) | `scripts/ComPlEx/ones_and_zeros.R` |
| `ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS`, `conserved_genes_COMPLETE_wider.RDS` (every complete clique of every orthogroup; no orthogroup excluded, no cap) | `scripts/ComPlEx/complete_cliques.R` |
| `results/orthogroup_size_rule.tsv`, `_model.tsv` (the even-odds size rule, D026: per eligible orthogroup, expressed genes, complete cliques, the fitted probability of a complete clique given size, and whether it is counted as conserved) | `scripts/ComPlEx/size_rule.R` (optional second argument: a fixed size threshold, for the sensitivity runs) |
| `ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS` (the conserved core: one complete clique per orthogroup counted by the size rule) and `conserved_genes_COMPLETE_filtered_flagged.RDS` (the orthogroups the rule does not count) | `scripts/ComPlEx/core_cliques.R <out> results/orthogroup_size_rule.tsv` |
| `ComPlEx/cliques/clique_HMcluster_list.rda` (the four clique clusters, analysis numbering) | `scripts/ComPlEx/clique_clusters.R` (with `clique_matrix.R`, shared with `make_Fig5.R`) |
| `DEGs/og_summary.RData` (per orthogroup, the genotypes with a gene and with a DEG) | `scripts/DE/og_summary.R` |

Run them in this order: `ones_and_zeros.R` and `complete_cliques.R`; `size_rule.R`; `core_cliques.R`, then
`clique_clusters.R`, before `make_TableS4.R`, `make_Fig5.R`, `make_Fig6.R` and `make_core_clusters.R`; `og_summary.R`
before `make_Fig1.R`, `make_TableS1.R` and `common_orthogroups.R`. Each writes its object to the path above, or to
the path given as its first argument. `scripts/tests/test_rebuilt_objects.R` runs the four into a temporary folder
and checks that each result is identical to the copy under `data/`.

`scripts/tests/check_hardcoded_results.py <value table>` searches the code for result values written as literals
(a result must be read from its producer's output). Its value table is the output of `scripts/analysis/value_changes.R`;
constants inspected and found not to be results are listed, with reasons, in `scripts/tests/hardcoded_allow.tsv`.

The size of each core orthogroup (expressed genes), its number of complete cliques, and the probability that an
orthogroup of its size is conserved in all 15 comparisons come from `scripts/validation/orthogroup_size_flag.R`
(`results/orthogroup_size_flag.tsv`, with the fitted model in `results/orthogroup_size_model.tsv`).

The orthology class of each orthogroup with a complete clique (D043) comes from the N1 (seed-plant) hierarchical
orthogroups inferred from multiple-alignment gene trees: `scripts/analysis/hog_clique_consistency.R <N1 table> <prefix>
<raw Orthogroups.tsv>` (per orthogroup and per representative clique), then `scripts/analysis/orthology_class.R <prefix>
<tree dir>` (`results/orthology_class.tsv`: orthologue- or paralogue-carried, the representative's consistency, the best
consistent clique, and which trees the call rests on). `scripts/analysis/orthology_class_tables.R` writes the descriptive
tables (class by size band, by manuscript cluster and by response pattern; paralogue-carried orthogroups with at most 20
genes). The class is a gene-tree call with error in both directions. `scripts/analysis/core_sensitivity.R` reports core
size, cluster sizes and TF cliques for working trees built with different size rules.

The per-species GO backgrounds (`annotation/bg_*.rda`) are the expressed-gene universe for each
species — every gene with nonzero VST expression across the cold time-course samples.
`scripts/annotation/make_backgrounds.R` regenerates them from the deposited expression tables.


Figure S6, panel a, comes from `scripts/validation/conservation_label_permutation.R` (deposit only).
Panel b needs the STRING v12 *A. thaliana* download (`3702.protein.links.detailed.v12.0.txt.gz` and
`3702.protein.aliases.v12.0.txt.gz` from string-db.org): run `make_core_clusters.R`, then
`string_external_validation.R`, then `make_FigS6.R` plots both. `string_by_cluster.R` (same inputs)
breaks the panel b test down by cluster and repeats it without manuscript cluster 2
(`results/string_by_cluster.tsv`, `results/string_within_cluster_pairs.tsv`).

Figure 1 (panel a, the UpSet of expressed-gene orthogroup overlap) is regenerated by
`make_Fig1.R` from two intermediate inputs provided in this repository, because they are not part
of the figshare deposit: the per-timepoint **median** expression tables
(`expression/med_expression/<sp>_med.txt.gz`) and the per-species gene→orthogroup maps
(`ComPlEx/gene2OG/gene2OG_<Sp>.rda`). A gene counts as expressed in a condition when its median
expression exceeds zero (control timepoint, or any cold timepoint); the shared expressed core is
the intersection across all twelve species×condition sets. Panels b and d (composition of expressed
genes and of DEGs in the shared orthogroups) reproduce exactly. Panel c (the DEG-category stacked
bars) classifies each DEG by how many species carry a gene, and a DEG, in its orthogroup — mapping
each DEG through the same gene→orthogroup maps and using the per-orthogroup `og_summary`
(categories: DEG across all six species / in some / species-specific / singletons) — and reproduces
the published counts exactly.

