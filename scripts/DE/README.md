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
