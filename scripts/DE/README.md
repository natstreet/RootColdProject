# Differential expression

DEGs in this study are defined as **adjusted P ≤ 0.05 and |log2 fold change| ≥ 1 (a two-fold change)**.

Method (per species): a `DESeqDataSet` is built from the raw gene-level count matrix with the
design `~ Condition`, where `Condition` is the single pooled control plus the four cold
timepoints (6 h, 24 h, 3 d, 10 d at 5 °C). After `DESeq()`, each cold timepoint is contrasted
against the control with

```r
results(dds, contrast = c("Condition", <timepoint>, "control"), alpha = 0.05, lfcThreshold = 0)
```

and genes with `padj ≤ 0.05` **and** `abs(log2FoldChange) >= 1` are differentially expressed
(direction from the sign of `log2FoldChange`). The per-species DEG set is the union over the four
timepoints (Col-0 5689, Ost-0 7995, aspen 11729, birch 3312, spruce 6353, pine 4717); these totals
are used for Figure 1 and Figure 2.

## Deposited files

Use these for the DEGs the paper reports:

- `DEGs/DEGs_<species>.RData` — per-species union DEG sets.
- `DEGs/per_timepoint_DEG_lists/<species>/…_p0.05_lfc0.txt` — the per-timepoint DEG lists.
- `DEGs/log2FoldChange_<code>.RData` — per-timepoint log2 fold changes, used for up/down direction.

`differential_expression.R` regenerates the per-timepoint lists and counts from the deposited
`data/dds/dds_<species>.rda` objects.
