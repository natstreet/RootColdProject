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

Outputs are written to `figures/` or `tables/`. R package requirements: `DESeq2`,
`matrixStats`, `topGO`, `GO.db`, `ggplot2`, `corrplot`, `dplyr`.

## Definitions used throughout

- **DEG:** adjusted P ≤ 0.05 and |log2 fold change| ≥ 1 (a two-fold change). (Deposited
  per-timepoint lists: `..._lfc0.txt` = the standard Wald test, `results()` with
  `lfcThreshold = 0`, then a post-hoc |log2FC| ≥ 1 filter — these are the sets used
  throughout the paper (every gene in them has |log2FC| ≥ 1); `..._l2fc1.txt` = the
  stricter fold-change-thresholded test (`lfcThreshold = 1`), which is not used.)
- **GO enrichment:** topGO, **classic** Fisher's exact test, biological process,
  Benjamini–Hochberg adjustment, padj < 0.05 (Tables S1, S2 and S3 all use the
  classic algorithm).

## Producer map and clean-room status (2026-08-06)

Status is from a clean-room run: unpack the deposit into `data/`, run each script,
nothing else on the path. "runs-from-deposit" = the script completes and writes its
artifact using only files in the archive.

| Output | Script | Key inputs (under `data/`) | Status |
|---|---|---|---|
| Fig. 2 | `scripts/figures/make_Fig2.R` | `DEGs/per_timepoint_DEG_lists/` | runs-from-deposit |
| Fig. 3C | `scripts/figures/make_Fig3C.R` | `*/Orthogroups_20240613.tsv`, `superclusters/…/SC*_gene_list.csv`, `superclusters/SC_mapping.csv`, `TranscriptionFactors/` | runs-from-deposit |
| Fig. S1 (matrix) | `scripts/figures/make_FigS1.R` | `superclusters/merged_super_clusters_means.RData` | runs-from-deposit¹ |
| Table S1 | `scripts/GO/make_TableS1.R` | `DEGs/og_summary.RData`, `DEGs/DEGs_col0.RData`, `*/Orthogroups_20240613.tsv`, `annotation/Athal_go_ids.tsv`, `annotation/bg_col0_ATC.rda` | runs-from-deposit |
| Table S2 | `scripts/GO/make_TableS2.R` | `ComPlEx/co_expressologs.RData`, `annotation/{go_ids,bg_*}`, `DEGs/DEGs_*.RData` | runs-from-deposit |
| Table S3 | `scripts/GO/make_TableS3.R` | `ComPlEx/cliques/{clique_HMcluster_list.rda,clique_genes_filterable_COMPLETE.RDS}`, `annotation/{go_ids,bg_*}` | runs-from-deposit² |
| Fig. 1 / Fig. 2 DEGs | `scripts/DE/differential_expression.R` | `data/dds/dds_{col0_soil,ost0_soil,Pt_new,Bp_new,Pa,Ps}.rda` | runs-from-deposit³ |

¹ Fig. S1: the script reproduces the 36×36 correlation matrix. The twelve
cross-species superclusters drawn as boxes on the published figure are a **manual**
grouping (Pearson r ≥ 0.7) and are not recoverable algorithmically (neither an
hclust cutree nor gene-list overlap reproduces it); their memberships are given by
`superclusters/SC_mapping.csv` and the per-SC gene lists.

² Table S3: regenerates the deposited table to Jaccard 0.963 (per-cell GO-set
overlap; shared-term p-values match to median Δlog10 0.011; residual ±1–2 terms per
cell is GO.db annotation-version drift). The spreadsheet p-value column
(`classicFisher`) holds topGO classic-Fisher values.

³ Differential expression: runs from the deposited raw `DESeqDataSet` objects
(`data/dds/`, integer counts + `~Condition` design). It reproduces the per-species
union DEG totals to within ~1.5–2.5% (col0 5544 vs 5689, ost0 7873 vs 7995, aspen
11534 vs 11729, birch 3277 vs 3312, spruce 6277 vs 6353, pine 4600 vs 4717); the
residual is the DESeq2 version used for independent filtering/dispersion (the
deposited lists were built with the 2025 release — record `sessionInfo()`). The
published per-timepoint and union DEG lists are also deposited
(`DEGs/DEGs_*.RData`, `DEGs/per_timepoint_DEG_lists/`) as the exact reference.

## Deposit consistency check

`scripts/DE/check_deposit_consistency.R` (run from a directory containing the unpacked deposit
`data/`) asserts the DEG-list invariants and prints a PASS/FAIL table: every `_lfc0`/`_lfc0585`/
`_l2fc1` list is |log2FC| ≥ 1; the per-species `_lfc0` unions equal `DEGs_*.RData`; thresholding
`padj_*.RData` at ≤ 0.05 reproduces the `_l2fc1` lists (so `padj_*` is the lfcThreshold = 1 test,
not the `_lfc0` DEGs); and `og_summary` all-six sums to 270. See `DEPOSIT_CHANGES.md` and
`scripts/DE/README.md` for what each deposited DEG artefact is.
