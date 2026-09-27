#!/usr/bin/env Rscript
# genotype_colours.R — the one genotype palette of the article, its figures, the deposit's explorer and exCompare.
#
# Colour-blind safe (Okabe–Ito based) and paired by lineage: the two Arabidopsis accessions in blues, the two broadleaf
# trees in greens, the two conifers in orange and pink. Every producer that colours by genotype sources this file and
# defines no genotype colours of its own (scripts/figures/check_genotype_colours.R enforces that).
#
#   source("scripts/figures/genotype_colours.R")   # defines GENOTYPE_COLOURS and genotype_tint()
#   Rscript scripts/figures/genotype_colours.R     # also writes scripts/figures/genotype_colours.json for the browser

GENOTYPE_COLOURS <- c(col0 = "#0072B2", ost0 = "#56B4E9", aspen = "#009E73", birch = "#7FBF3F",
                      spruce = "#D55E00", pine = "#CC79A7")

# a lighter version of each colour, mixed with white; for filled areas that carry black text
genotype_tint <- function(col, amount = 0.5) {
  m <- grDevices::col2rgb(col) / 255
  m <- m + (1 - m) * amount
  stats::setNames(grDevices::rgb(m[1, ], m[2, ], m[3, ]), names(col))
}

if (sys.nframe() == 0) {
  out <- "scripts/figures/genotype_colours.json"
  writeLines(paste0("{\n", paste(sprintf('  "%s": "%s"', names(GENOTYPE_COLOURS), GENOTYPE_COLOURS), collapse = ",\n"),
                    "\n}"), out)
  cat("wrote", out, "\n")
}
