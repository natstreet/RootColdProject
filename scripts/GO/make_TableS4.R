#!/usr/bin/env Rscript
# make_TableS4.R — Table S4: gene membership of the four super clusters, per species/ecotype.
# One worksheet per super cluster (SC_I–SC_IV); within a worksheet, one column per
# species/ecotype listing the differentially expressed genes assigned to that super cluster in
# that species (genes whose cold-response profile correlates with the super cluster mean at
# Pearson r >= 0.7). A final worksheet, "TF_families", lists the transcription factors among the
# members (all species pooled in one column) with the family each belongs to.
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs:
#   data/superclusters/membership/<species>/cluster_<n>_members.csv   per-species SC membership
#   data/TranscriptionFactors/TF_<species>.RData                      gene -> TF family
# Output: Table_S4_supercluster_gene_membership.xlsx
#
# Usage: Rscript scripts/GO/make_TableS4.R

suppressPackageStartupMessages({ library(openxlsx) })

SC <- c(SC_I = 6, SC_II = 2, SC_III = 3, SC_IV = 1)          # main-text label -> original number
species <- c(col0 = "Col-0", ost0 = "Ost-0", birch = "Birch",
             aspen = "Aspen", pine = "Pine", spruce = "Spruce")

ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
famtab <- function(sp){
  x <- ld(sprintf("data/TranscriptionFactors/TF_%s.RData", if (sp == "ost0") "col0" else sp))
  if (all(c("Gene_ID","Family") %in% colnames(x)))
    setNames(as.character(x$Family), as.character(x$Gene_ID))
  else setNames(as.character(x[[2]]), as.character(x[[1]]))
}
fam <- unlist(lapply(names(species), famtab)); fam <- fam[!duplicated(names(fam))]

members <- function(sp, n)
  sort(read.csv(sprintf("data/superclusters/membership/%s/cluster_%d_members.csv", sp, n),
                stringsAsFactors = FALSE)$gene_id)

wb  <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom")
tf_seen <- character(0)
for (lab in names(SC)) {
  cols <- lapply(names(species), members, n = SC[[lab]]); names(cols) <- unname(species)
  n <- max(lengths(cols), 1L)
  df <- as.data.frame(lapply(cols, function(x) c(x, rep(NA_character_, n - length(x)))),
                      check.names = FALSE, stringsAsFactors = FALSE)
  addWorksheet(wb, lab); writeData(wb, lab, df, headerStyle = hdr, keepNA = FALSE)
  setColWidths(wb, lab, cols = seq_along(df), widths = 16); freezePane(wb, lab, firstRow = TRUE)
  tf_seen <- union(tf_seen, unlist(cols)[unlist(cols) %in% names(fam)])
}
# TF_families worksheet: every TF among the members, with its family
tf <- data.frame(Gene = sort(tf_seen), `TF family` = fam[sort(tf_seen)],
                 check.names = FALSE, row.names = NULL)
addWorksheet(wb, "TF_families"); writeData(wb, "TF_families", tf, headerStyle = hdr)
setColWidths(wb, "TF_families", cols = 1:2, widths = c(16, 14)); freezePane(wb, "TF_families", firstRow = TRUE)

saveWorkbook(wb, "Table_S4_supercluster_gene_membership.xlsx", overwrite = TRUE)

for (lab in names(SC))
  cat(sprintf("%-6s: %s\n", lab,
      paste(sprintf("%s=%d", unname(species),
                    vapply(names(species), function(sp) length(members(sp, SC[[lab]])), integer(1))),
            collapse = "  ")))
cat(sprintf("TF_families: %d transcription factors across %d families\n",
            nrow(tf), length(unique(tf$`TF family`))))
