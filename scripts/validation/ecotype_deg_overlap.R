#!/usr/bin/env Rscript
# ecotype_deg_overlap.R — how many Col-0 DEGs are also Ost-0 DEGs.
#
# Each ecotype's DEG set is the union over the four timepoints of its per-timepoint DEG lists (the *_p0.05_lfc0.txt
# files: padj <= 0.05 against a null of no change and |log2 fold change| >= 1). The union must equal the published set
# in data/DEGs/DEGs_<ecotype>.RData; the script stops otherwise.
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/DEGs/per_timepoint_DEG_lists/{col0,ost0}/DEGs_<ecotype>_C_vs_<tp>_p0.05_lfc0.txt
#   data/DEGs/DEGs_{col0,ost0}.RData
# Output: results/ecotype_deg_overlap.tsv (quantity, value)
#
# Usage: Rscript scripts/validation/ecotype_deg_overlap.R

ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
degs <- function(eco) {
  f <- Sys.glob(sprintf("data/DEGs/per_timepoint_DEG_lists/%s/DEGs_%s_C_vs_*_p0.05_lfc0.txt", eco, eco))
  stopifnot(length(f) == 4)
  u <- unique(unlist(lapply(f, function(x) { g <- trimws(readLines(x)); g[nzchar(g)] })))
  pub <- unique(as.character(ld(sprintf("data/DEGs/DEGs_%s.RData", eco))))
  if (!setequal(u, pub)) stop(eco, ": union of the per-timepoint lists differs from DEGs_", eco, ".RData")
  u
}
col0 <- degs("col0"); ost0 <- degs("ost0"); both <- intersect(col0, ost0)
out <- data.frame(quantity = c("col0_DEGs", "ost0_DEGs", "col0_DEGs_also_ost0_DEGs", "pct_of_col0_DEGs",
                               "ost0_DEGs_also_col0_DEGs_pct"),
                  value = c(length(col0), length(ost0), length(both), 100 * length(both) / length(col0),
                            100 * length(both) / length(ost0)))
dir.create("results", showWarnings = FALSE)
write.table(out, "results/ecotype_deg_overlap.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
