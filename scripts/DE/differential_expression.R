#!/usr/bin/env Rscript
# differential_expression.R — per-species differential expression from the count matrix.
#
# DEG definition: adjusted P <= 0.05 AND |log2 fold change| >= 1 (a two-fold change).
#
# Method (per species): starting from the raw gene-level DESeqDataSet, DESeq2's negative-binomial
# Wald test contrasts each cold timepoint (6h, 24h, 3d, 10d at 5 C) against the pooled control,
# with the median normalised count as the independent-filter statistic (UPSCb `extract_results`
# convention). Genes with padj <= 0.05 and |log2FoldChange| >= 1 are differentially expressed;
# direction is the sign of log2FoldChange. Per-timepoint unions give the per-species DEG totals
# (Col-0 5689, Ost-0 7995, aspen 11729, birch 3312, spruce 6353, pine 4717).
#
# NOTE on the deposited list names: `..._lfc0.txt` are the DEG sets used throughout the paper —
# the standard Wald test (results() with lfcThreshold = 0) with a post-hoc |log2FC| >= 1 filter
# (every gene in them has |log2FC| >= 1). `..._l2fc1.txt` instead use a fold-change-thresholded
# null (results() with lfcThreshold = 1) and are the stricter, unused alternative. The gene sets
# reproduce to Jaccard ~0.96; a residual ~3-7% is the DESeq2-version dependence of
# dispersion/independent-filtering (deposited built with the 2025 release; record sessionInfo()).
#
# Inputs (unpack data.tar.gz from figshare into ./data; the dds objects must be added there):
#   data/dds/dds_<species>.rda   (raw DESeqDataSet per species)
# Output:
#   DE/<species>/DEGs_<species>_C_vs_<tp>.txt  and  DE/<species>/DEG_counts.tsv
#
# Usage: Rscript scripts/DE/differential_expression.R

suppressMessages({ library(DESeq2); library(matrixStats) })

# species -> (dds file, control-level label). Design variable is read from each dds.
SP <- list(
  col0   = c(dds = "dds_col0_soil", ctrl = "control"),
  ost0   = c(dds = "dds_ost0_soil", ctrl = "control"),
  aspen  = c(dds = "dds_Pt_new",    ctrl = "control"),
  birch  = c(dds = "dds_Bp_new",    ctrl = "control"),
  spruce = c(dds = "dds_Pa",        ctrl = "control"),
  pine   = c(dds = "dds_Ps",        ctrl = "control"))

getdds <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

for (sp in names(SP)) {
  f <- Sys.glob(sprintf("data/dds/%s.rda", SP[[sp]]["dds"]))
  if (!length(f)) { message("skip ", sp, " (dds not found)"); next }
  dds <- DESeq(getdds(f[1]), quiet = TRUE)
  dv  <- all.vars(design(dds))[1]
  ctrl <- SP[[sp]]["ctrl"]
  levs <- setdiff(levels(colData(dds)[[dv]]), ctrl)
  mf  <- rowMedians(counts(dds, normalized = TRUE))
  outdir <- file.path("DE", sp); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  counts_tab <- data.frame(); u <- character(0)
  for (lv in levs) {
    r <- results(dds, contrast = c(dv, lv, ctrl), alpha = 0.05, lfcThreshold = 0, filter = mf)
    deg <- rownames(r)[!is.na(r$padj) & r$padj <= 0.05 & abs(r$log2FoldChange) >= 1]
    up  <- sum(r[deg, "log2FoldChange"] > 0); down <- length(deg) - up
    writeLines(deg, file.path(outdir, sprintf("DEGs_%s_C_vs_%s.txt", sp, lv)))
    counts_tab <- rbind(counts_tab, data.frame(timepoint = lv, up = up, down = down, total = length(deg)))
    u <- union(u, deg)
  }
  write.table(counts_tab, file.path(outdir, "DEG_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%-7s union DEGs = %d\n", sp, length(u)))
}
