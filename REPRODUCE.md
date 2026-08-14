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
`matrixStats`, `topGO`, `GO.db`, `ggplot2`, `corrplot`, `dplyr`, `tidyr`, `openxlsx`.

## Definitions used throughout

- **DEG:** adjusted P ≤ 0.05 and |log2 fold change| ≥ 1 (a two-fold change). The deposited
  per-timepoint lists (`…_p0.05_lfc0.txt`) and per-species union sets (`DEGs_*.RData`) are the
  DEGs used throughout the paper.
- **GO enrichment:** topGO, classic Fisher's exact test, biological process,
  Benjamini–Hochberg adjustment, padj < 0.05 (Tables S1, S2 and S3).
- **Super clusters:** main-text labels SC_I–SC_IV map to the original super cluster numbers
  SC_I = 6, SC_II = 2, SC_III = 3, SC_IV = 1 (`superclusters/SC_mapping.csv`).

## Producer map

| Output | Script | Key inputs (under `data/`) |
|---|---|---|
| Fig. 1 | `scripts/figures/make_Fig1.R` | `dds/`, `expression/`, `ComPlEx/Orthogroups_20240613.tsv`, `DEGs/` |
| Fig. 2 | `scripts/figures/make_Fig2.R` | `DEGs/per_timepoint_DEG_lists/`, `DEGs/log2FoldChange_*.RData` |
| Fig. 3C | `scripts/figures/make_Fig3C.R` | `*/Orthogroups_20240613.tsv`, `superclusters/…/SC*_gene_list.csv`, `superclusters/SC_mapping.csv`, `TranscriptionFactors/` |
| Fig. S1 (matrix) | `scripts/figures/make_FigS1.R` | `superclusters/merged_super_clusters_means.RData` |
| Fig. S2 | `scripts/figures/make_FigS2.R` | `superclusters/tf_membership/`, `TranscriptionFactors/` |
| Table S1 | `scripts/GO/make_TableS1.R` | `DEGs/og_summary.RData`, `DEGs/DEGs_col0.RData`, `*/Orthogroups_20240613.tsv`, `annotation/Athal_go_ids.tsv`, `annotation/bg_col0_ATC.rda` |
| Table S2 | `scripts/GO/make_TableS2.R` | `ComPlEx/co_expressologs.RData`, `annotation/{go_ids,bg_*}`, `DEGs/DEGs_*.RData` |
| Table S3 | `scripts/GO/make_TableS3.R` | `ComPlEx/cliques/{clique_HMcluster_list.rda,clique_genes_filterable_COMPLETE.RDS}`, `annotation/{go_ids,bg_*}` |
| Table S4 | `scripts/GO/make_TableS4.R` | `superclusters/tf_membership/` |
| DEGs (Fig. 1 / Fig. 2) | `scripts/DE/differential_expression.R` | `dds/dds_{col0_soil,ost0_soil,Pt_new,Bp_new,Pa,Ps}.rda` |

The twelve cross-species super clusters drawn on Figure S1 are a manual grouping of the
per-species clusters (Pearson r ≥ 0.7); their memberships are given by
`superclusters/SC_mapping.csv` and the per-SC gene lists. Figure S2 and Table S4 draw on the
per-species transcription-factor membership of each super cluster in
`superclusters/tf_membership/`.
