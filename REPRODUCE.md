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
   Rscript scripts/GO/make_TableS3.R
   ```

Outputs are written to the working directory. R package requirements: `DESeq2`,
`matrixStats`, `topGO`, `GO.db`, `ggplot2`, `ggtext`, `corrplot`, `dplyr`, `tidyr`,
`pheatmap`, `patchwork`, `UpSetR`, `readr`, `openxlsx`.

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
  Benjamini–Hochberg adjustment, padj < 0.05 (Tables S1, S2 and S3).
- **Super clusters:** main-text labels SC_I–SC_IV map to the original super cluster numbers
  SC_I = 6, SC_II = 2, SC_III = 3, SC_IV = 1 (`superclusters/SC_mapping.csv`).

## Producer map

| Output | Script | Key inputs (under `data/`) |
|---|---|---|
| Fig. 1 | `scripts/figures/make_Fig1.R` | `expression/med_expression/`, `ComPlEx/gene2OG/`, `annotation/bg_*`, `DEGs/`; panel values in `results/Figure1_panel_values.tsv` |
| Fig. 2 | `scripts/figures/make_Fig2.R` | `DEGs/per_timepoint_DEG_lists/`, `DEGs/log2FoldChange_*.RData` |
| Fig. 3 | `scripts/figures/make_Fig3.R` | `superclusters/membership/`, `expression/`, `superclusters/DEGs_TFs/`, `*/Orthogroups_20240613.tsv`, `superclusters/SC_mapping.csv`, `TranscriptionFactors/` |
| Fig. 4 | `scripts/figures/make_Fig4.R` | `ComPlEx/co_expressologs.RData`, `ComPlEx/all_go_results_clades.RData`; Venn regions in `results/Figure4_venn_regions.tsv` |
| Fig. 5 | `scripts/figures/make_Fig5.R` | `ComPlEx/cliques/{conserved_genes_COMPLETE_filtered.RDS,clique_HMcluster_list.rda}`, `expression/med_expression/`; the Table S3 output and `scripts/figures/Fig5_captions.tsv` (block captions, by GO ID) |
| Fig. 6 | `scripts/figures/make_Fig6.R` | `annotation/gene_aliases_20140331.txt`, `ComPlEx/cliques/{clique_genes_filterable_COMPLETE.RDS,conserved_genes_COMPLETE_filtered.RDS}`, `expression/`, `DEGs/` |
| Fig. S1 | `scripts/figures/make_FigS1.R` | `superclusters/merged_super_clusters_means.RData`, `superclusters/SC_mapping.csv` |
| Fig. S2 | `scripts/figures/make_FigS2.R` | `superclusters/membership/`, `TranscriptionFactors/` |
| Fig. S3 | `scripts/figures/make_FigS3.R` | `ComPlEx/co_expressologs.RData`, `DEGs/DEGs_*.RData` |
| Fig. S4 | `scripts/figures/make_FigS4.R` | `ComPlEx/co_expressologs.RData` |
| Fig. S5 | `scripts/figures/make_FigS5.py` | outputs of `scripts/validation/conservation_label_permutation.R` (panel a) and `make_core_clusters.R` + `string_external_validation.py` (panel b; **STRING download, see note**) |
| Fig. S6 | `scripts/figures/make_FigS6.R` | `ComPlEx/co_expressologs.RData`, `expression/` |
| Fig. S7 | `scripts/figures/make_FigS7.py` | `open/{conserved_core,gene_index,co_expressologs,transcription_factors}.tsv`, `networks/network_*_hrr100.tsv` — **deposit v4 or later**; `DE/<species>/` from `differential_expression.R` |
| Table S1 (`Table_S1.xlsx`) | `scripts/GO/make_TableS1.R` (with `go_table_utils.R`) | `DEGs/og_summary.RData`, `DEGs/DEGs_col0.RData`, `*/Orthogroups_20240613.tsv`, `annotation/Athal_go_ids.tsv`, `annotation/bg_col0_ATC.rda` |
| Table S2 (`Table_S2.xlsx`) | `scripts/GO/make_TableS2.R` (with `go_table_utils.R`) | `ComPlEx/co_expressologs.RData`, `annotation/{go_ids,bg_*}`, `DEGs/DEGs_*.RData` |
| Table S3 (`Table_S3.xlsx`) | `scripts/GO/make_TableS3.R` (with `go_table_utils.R`) | `ComPlEx/cliques/{clique_HMcluster_list.rda,clique_genes_filterable_COMPLETE.RDS}`, `annotation/{go_ids,bg_*}` |
| Table S4 | `scripts/GO/make_TableS4.R` | `superclusters/membership/` |
| DEGs (Fig. 1 / Fig. 2) | `scripts/DE/differential_expression.R` | `dds/dds_{col0_soil,ost0_soil,Pt_new,Bp_new,Pa,Ps}.rda` |
| GO backgrounds | `scripts/annotation/make_backgrounds.R` | `expression/<sp>_root_cold_stress_expression.txt` |
| Species networks | `scripts/networks/build_species_networks.R` | `expression/<sp>_root_cold_stress_expression.txt` |
| TF clique responses (Discussion) | `scripts/analysis/tf_clique_directions.R` | `dds/`, `expression/`, `ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS`, `annotation/gene_aliases_20140331.txt`; output `results/tf_clique_directions.tsv` |
| Orthogroups split into several HOGs (Methods) | `scripts/analysis/orthogroup_hog_split.R` | `ComPlEx/Orthogroups_20240613.tsv`, `ComPlEx/Phylogenetic_Hierarchical_Orthogroups_N0.tsv.gz` (in this repository; N0.tsv of the same OrthoFinder run, from github.com/tallgran/conifer-JasonHill); output `results/orthogroup_hog_split.tsv` |
| Cluster profiles (Results, clusters 1–2) | `scripts/analysis/cluster_profiles.py` | `open/{gene_index,conserved_core,orthogroup_annotation}.tsv` — **deposit v4 or later**; outputs in `results/cluster_profiles/` |

The per-species GO backgrounds (`annotation/bg_*.rda`) are the expressed-gene universe for each
species — every gene with nonzero VST expression across the cold time-course samples.
`scripts/annotation/make_backgrounds.R` regenerates them from the deposited expression tables.

Figure S7 reads the open-format tables in `data/open/` and the per-species co-expression networks
in `data/networks/`, both introduced in **version 4** of the figshare deposit. It will not run
against versions 1–3.

Figure S5, panel a, comes from `scripts/validation/conservation_label_permutation.R` (deposit only).
Panel b needs the STRING v12 *A. thaliana* download (`3702.protein.links.detailed.v12.0.txt.gz` and
`3702.protein.aliases.v12.0.txt.gz` from string-db.org): run `make_core_clusters.R`, then
`string_external_validation.py`, then `make_FigS5.py` plots both. `string_by_cluster.py` (same inputs)
breaks the panel b test down by cluster and repeats it without manuscript cluster 2
(`results/string_by_cluster.tsv`, `results/string_within_cluster_pairs.tsv`). (`conservation_null_test.py` is a
separate network-rewiring test and is not what Figure S5 shows.)

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
files, and `superclusters/SC_mapping.csv` names the four discussed in the main text. Figure S2 and Table S4 draw on the
per-species, per-ecotype gene membership of each super cluster in `superclusters/membership/`
(the differentially expressed genes assigned to a super cluster in a species, i.e. those whose
cold-response profile correlates with the super cluster mean at r ≥ 0.7); this membership is
included in the repository and is also provided as Supplementary Table S4.
