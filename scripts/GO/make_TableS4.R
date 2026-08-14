#!/usr/bin/env Rscript
# make_TableS4.R — Table S4: per-species transcription-factor membership of the four super
# clusters. One worksheet per super cluster (SC_I–SC_IV); within a worksheet, one column per
# species/ecotype listing the transcription-factor genes assigned to that super cluster in that
# species (the TF-DEC genes summarised by Figure S2).
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Input:  data/superclusters/tf_membership/<species>/cluster_<n>_TFDEC_members.csv
# Output: Table_S4_supercluster_TF_membership.xlsx
#
# Usage: Rscript scripts/GO/make_TableS4.R

suppressPackageStartupMessages({ library(openxlsx) })

SC <- c(SC_I = 6, SC_II = 2, SC_III = 3, SC_IV = 1)          # main-text label -> original number
species <- c(col0 = "Col-0", ost0 = "Ost-0", birch = "Birch",
             aspen = "Aspen", pine = "Pine", spruce = "Spruce")

members <- function(sp, n){
  f <- sprintf("data/superclusters/tf_membership/%s/cluster_%d_TFDEC_members.csv", sp, n)
  sort(read.csv(f, stringsAsFactors = FALSE)$gene_id)
}

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", border = "bottom")
for (lab in names(SC)) {
  cols <- lapply(names(species), members, n = SC[[lab]])
  names(cols) <- unname(species)
  n <- max(lengths(cols), 1L)
  df <- as.data.frame(lapply(cols, function(x) c(x, rep(NA_character_, n - length(x)))),
                      check.names = FALSE, stringsAsFactors = FALSE)
  addWorksheet(wb, lab)
  writeData(wb, lab, df, headerStyle = hdr, keepNA = FALSE)
  setColWidths(wb, lab, cols = seq_along(df), widths = 16)
  freezePane(wb, lab, firstRow = TRUE)
}
saveWorkbook(wb, "Table_S4_supercluster_TF_membership.xlsx", overwrite = TRUE)

for (lab in names(SC))
  cat(sprintf("%-6s: %s\n", lab,
      paste(sprintf("%s=%d", unname(species),
                    vapply(names(species), function(sp) length(members(sp, SC[[lab]])), integer(1))),
            collapse = "  ")))
