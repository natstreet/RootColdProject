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
                      spruce = "#D55E00", pine = "#CC79A7",
                      # two more for the other studies in PlantGenIE's comparative layer (Rodriguez et al. 2026: wild
                      # cherry, a broadleaf, in the green group; lodgepole pine, a conifer, in the warm group). The
                      # article's figures do not use them. Chosen so that each is further from every other colour
                      # (CIELAB difference at least 17.3, under normal vision and simulated deutan, protan and tritan
                      # vision) than the closest pair of the original six (16.1, Ost-0 and aspen under tritan vision).
                      # The Rodriguez et al. figures colour lineages, not species, so there was no species colour to reuse.
                      cherry = "#4A6122", lodgepole = "#E62E5C")

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
