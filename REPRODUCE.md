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
(`DEGs/og_summary.RData`), the per-timepoint DEG lists (`DEGs/per_timepoint_DEG_lists/`), the
super-cluster mapping and mean profiles (`superclusters/SC_mapping.csv`,
`superclusters/merged_super_clusters_means.RData`) and the super-cluster memberships
(`superclusters/membership/`). Unpacking the deposit alongside them gives a complete `data/` tree.

## Definitions used throughout

- **DEG:** adjusted P ≤ 0.05 and |log2 fold change| ≥ 1 (a two-fold change). The per-timepoint
  lists (`…_p0.05_lfc0.txt`, provided in this repository) and the deposited per-species union sets
  (`DEGs_*.RData`) are the DEGs used throughout the paper.
- **GO enrichment:** topGO, classic Fisher's exact test, biological process,
  Benjamini–Hochberg adjustment, padj < 0.05 (Tables S1, S3 and S4).
- **Super clusters:** main-text labels SC_I–SC_IV map to the original super cluster numbers
  SC_I = 6, SC_II = 2, SC_III = 3, SC_IV = 1 (`superclusters/SC_mapping.csv`).

## Producer map

| Output | Script | Key inputs (under `data/`) |
|---|---|---|
| Fig. 1 | `scripts/figures/make_Fig1.R` | `expression/med_expression/`, `ComPlEx/gene2OG/`, `annotation/bg_*`, `DEGs/`; panel values in `results/Figure1_panel_values.tsv` |
| Fig. 2 | `scripts/figures/make_Fig2.R` | `DEGs/per_timepoint_DEG_lists/`, `DEGs/log2FoldChange_*.RData` |
| Fig. 3 | `scripts/figures/make_Fig3.R` | `superclusters/membership/`, `expression/`, `superclusters/DEGs_TFs/`, `*/Orthogroups_20240613.tsv`, `superclusters/SC_mapping.csv`, `TranscriptionFactors/` |
| Fig. 4 | `scripts/figures/make_Fig4.R` | `ComPlEx/co_expressologs.RData`; panel b from the Table S3 output (`tables/TableS3_cladeGO_classic_p10.tsv`); Venn regions in `results/Figure4_venn_regions.tsv`, panel b values in `results/Figure4b_values.tsv` |
| Fig. 5 | `scripts/figures/make_Fig5.R` | `ComPlEx/cliques/{conserved_genes_COMPLETE_filtered.RDS,clique_HMcluster_list.rda}`, `expression/med_expression/`; the Table S4 output and `scripts/figures/Fig5_captions.tsv` (block captions, by GO ID) |
| Fig. 6 | `scripts/figures/make_Fig6.R` | `annotation/gene_aliases_20140331.txt`, `ComPlEx/cliques/{clique_genes_filterable_COMPLETE.RDS,conserved_genes_COMPLETE_filtered.RDS}`, `expression/`, `DEGs/` |
| Fig. S1 | `scripts/figures/make_FigS1.R` | `superclusters/merged_super_clusters_means.RData`, `superclusters/SC_mapping.csv` |
| Fig. S2 | `scripts/figures/make_FigS2.R` | `open/{conserved_core,gene_index,co_expressologs,transcription_factors}.tsv`, `networks/network_*_hrr100.tsv` — **deposit v4 or later**; `DE/<species>/` from `differential_expression.R` |
| Fig. S3 | `scripts/figures/make_FigS3.R` | `superclusters/membership/`, `TranscriptionFactors/` |
| Fig. S4 | `scripts/figures/make_FigS4.R` | `ComPlEx/co_expressologs.RData`, `DEGs/DEGs_*.RData` |
| Fig. S5 | `scripts/figures/make_FigS5.R` | `ComPlEx/co_expressologs.RData`; per-comparison counts in `results/Figure_S5_values.tsv` |
| Fig. S6 | `scripts/figures/make_FigS6.R` | outputs of `scripts/validation/conservation_label_permutation.R` (panel a) and `make_core_clusters.R` + `string_external_validation.R` (panel b; **STRING download, see note**) |
| Fig. S7 | `scripts/figures/make_FigS7.R` | `ComPlEx/co_expressologs.RData`, `expression/` |
| Table S1 (`Table_S1.xlsx`) | `scripts/GO/make_TableS1.R` (with `go_table_utils.R`) | `DEGs/og_summary.RData`, `DEGs/DEGs_col0.RData`, `*/Orthogroups_20240613.tsv`, `annotation/Athal_go_ids.tsv`, `annotation/bg_col0_ATC.rda` |
| Table S2 (`Table_S2.xlsx`) | `scripts/GO/make_TableS2.R` | `superclusters/membership/` |
| Table S3 (`Table_S3.xlsx`) | `scripts/GO/make_TableS3.R` (with `go_table_utils.R`) | `ComPlEx/co_expressologs.RData`, `annotation/{go_ids,bg_*}`, `DEGs/DEGs_*.RData` |
| Table S4 (`Table_S4.xlsx`) | `scripts/GO/make_TableS4.R` (with `go_table_utils.R`) | `ComPlEx/cliques/{clique_HMcluster_list.rda,clique_genes_filterable_COMPLETE.RDS}`, `annotation/{go_ids,bg_*}` |
| DEGs (Fig. 1 / Fig. 2) | `scripts/DE/differential_expression.R` | `dds/dds_{col0_soil,ost0_soil,Pt_new,Bp_new,Pa,Ps}.rda` |
| GO backgrounds | `scripts/annotation/make_backgrounds.R` | `expression/<sp>_root_cold_stress_expression.txt` |
| Species networks | `scripts/networks/build_species_networks.R` | `expression/<sp>_root_cold_stress_expression.txt` |
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
| Deposit `data/open/` additions: GO tables, manuscript cluster numbers, DEG calls, ENA runs | `scripts/deposit/{export_go_tables,manuscript_clusters,de_calls,sample_runs_columns}.R` | Tables S1/S3/S4 outputs, `DE/de_results.tsv`, `results/sample_ena_runs.tsv`; each edits the deposit's `open/` tables in place |

The per-species GO backgrounds (`annotation/bg_*.rda`) are the expressed-gene universe for each
species — every gene with nonzero VST expression across the cold time-course samples.
`scripts/annotation/make_backgrounds.R` regenerates them from the deposited expression tables.

Figure S2 reads the open-format tables in `data/open/` and the per-species co-expression networks
in `data/networks/`, both introduced in **version 4** of the figshare deposit. It will not run
against versions 1–3.

Figure S6, panel a, comes from `scripts/validation/conservation_label_permutation.R` (deposit only).
Panel b needs the STRING v12 *A. thaliana* download (`3702.protein.links.detailed.v12.0.txt.gz` and
`3702.protein.aliases.v12.0.txt.gz` from string-db.org): run `make_core_clusters.R`, then
`string_external_validation.R`, then `make_FigS6.R` plots both. `string_by_cluster.R` (same inputs)
breaks the panel b test down by cluster and repeats it without manuscript cluster 2
(`results/string_by_cluster.tsv`, `results/string_within_cluster_pairs.tsv`). (`conservation_null_test.py` is a
separate network-rewiring test and is not what Figure S6 shows.)

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

The twelve cross-species super clusters drawn on Figure S1 are the groups obtained by hierarchical
clustering (Ward's method, on 1 − Pearson r) of the 36 per-species cluster-mean profiles, with the
tree cut into twelve; they are numbered 1–12 in dendrogram order, the numbering of the membership
files, and `superclusters/SC_mapping.csv` names the four discussed in the main text. Figure S3 and Table S2 draw on the
per-species, per-ecotype gene membership of each super cluster in `superclusters/membership/`
(the differentially expressed genes assigned to a super cluster in a species, i.e. those whose
cold-response profile correlates with the super cluster mean at r ≥ 0.7); this membership is
included in the repository and is also provided as Supplementary Table S2.
