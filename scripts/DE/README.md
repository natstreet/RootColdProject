# Differential expression

DEGs in this study are defined as **adjusted P ≤ 0.05 and |log2 fold change| ≥ 1 (a two-fold change)**.

Method (per species): a `DESeqDataSet` is built from the raw gene-level count matrix with the
design `~ Condition`, where `Condition` is the single pooled control plus the four cold
timepoints (6 h, 24 h, 3 d, 10 d at 5 °C). After `DESeq()`, each cold timepoint is contrasted
against the control with

```r
results(dds, contrast = c("Condition", <timepoint>, "control"), alpha = 0.05, lfcThreshold = 0)
```

and genes with `padj ≤ 0.05` **and** `abs(log2FoldChange) >= 1` are taken as differentially expressed
(direction from the sign of `log2FoldChange`). The per-timepoint gene lists are provided in the
figshare deposit under `DEGs/per_timepoint_DEG_lists/<species>/…_p0.05_lfc0.txt` (these carry the |log2FC| ≥ 1 cutoff despite the `lfc0` name); their per-species
unions are the `DEGs_*` totals (Col-0 5689, Ost-0 7995, aspen 11729, birch 3312, spruce 6353,
pine 4717) used for Figure 1 and Figure 2. Lists at |log2FC| ≥ 1 (`l2fc1`) and ≥ 0.585 (`lfc0585`)
are also provided for reference.

The original per-species DE scripts are by A. Vergara / T. Aro (`DEGs_tuuli/` in the analysis
working directory). Exact reproduction from the raw count matrix requires the per-species
`DESeqDataSet` objects, since the deposited DE used the gene set carried in those objects; a naïve
rebuild from the full unfiltered count matrix over-counts. The dds/count objects can be added to
the deposit on request.

## The three per-timepoint list families (read this before reusing them)

The filename suffix names the `lfcThreshold` argument passed to DESeq2 `results()`, **not** the
fold-change cutoff. All three families additionally apply a post-hoc `|log2FC| ≥ 1` filter, so
**every gene in every list (all three families) has |log2FC| ≥ 1**:

| suffix     | `results()` null       | meaning                                                              |
|------------|------------------------|---------------------------------------------------------------------|
| `_lfc0`    | `lfcThreshold = 0`     | standard Wald test + post-hoc \|log2FC\|≥1 — **the DEGs used throughout the paper** |
| `_lfc0585` | `lfcThreshold = 0.585` | 1.5-fold-thresholded null + post-hoc cut (stricter; reference only)  |
| `_l2fc1`   | `lfcThreshold = 1`     | two-fold-thresholded null (strictest; reference only)               |

`_lfc0` is what Figures 1, 2 and the `DEGs_*.RData` objects use. The name "lfc0" looks like
"no fold filter"; it is not.

## `padj_*.RData` — do NOT use to recheck the paper's DEGs or figure significance marks

`data/DEGs/padj_<code>.RData` holds adjusted P-values from the **`lfcThreshold = 1`** test.
Thresholding them at `padj ≤ 0.05` reproduces the strict `_l2fc1` lists (verified 24/24
species × timepoint), **not** the `_lfc0` DEGs the paper reports. A reuser checking whether a gene
shown in Figure 6 is "significant" against `padj_*.RData` gets the stricter answer and may wrongly
conclude the figure is off. To recheck a paper DEG use the `_lfc0` lists / `DEGs_*.RData`, and
`log2FoldChange_<code>.RData` for direction.

## `log2FoldChange_ps_n.RData` — not the published pine object

`data/DEGs/log2FoldChange_ps_n.RData` is a second pine log2FC object, same dimensions and gene
identifiers as `log2FoldChange_ps.RData` but substantively different values (≈40k/49k rows differ
at 6 h; correlation ≈0.15; the deposited pine DEG lists agree with `_ps`, not `_ps_n`). Its exact
provenance could not be established from the deposited files. **The published pine results use
`log2FoldChange_ps.RData`; `_ps_n` should not be used** and is a candidate for removal in a future
deposit version.

## Consistency check

`scripts/DE/check_deposit_consistency.R`, run from a directory containing the unpacked deposit
`data/`, asserts all of the above (every list is |log2FC|≥1; `_lfc0` unions equal `DEGs_*`;
`padj_*`@0.05 equals `_l2fc1`; `og_summary` all-six = 270) and prints a PASS/FAIL table.
