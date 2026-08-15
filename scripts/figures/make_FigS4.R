#!/usr/bin/env Rscript
# make_FigS4.R — Figure S4: UpSet plot of orthogroup overlaps between the co-expressolog sets of
# the pairwise species comparisons (MaxpVal < 0.1). Horizontal bars give the number of
# co-expressolog orthogroups per comparison; the matrix/vertical bars give the intersections.
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs: data/ComPlEx/co_expressologs.RData
# Output: Figure_S4.pdf
# Usage:  Rscript scripts/figures/make_FigS4.R

suppressPackageStartupMessages({ library(UpSetR); library(dplyr) })
P <- 0.1
ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
co <- ld("data/ComPlEx/co_expressologs.RData")

co$pair <- apply(cbind(co$Species1, co$Species2), 1, function(v) paste(sort(v), collapse = "_"))
d <- co %>% filter(MaxpVal < P) %>% distinct(OrthologGroup, pair)
M <- (as.data.frame.matrix(table(d$OrthologGroup, d$pair) > 0)) * 1
allp <- colnames(M)

teal <- "#2A9BAF"; orange <- "#C46E00"
q <- list(
  list(query = intersects, params = list("col0_ost0"),   color = teal,   active = TRUE),
  list(query = intersects, params = list("pine_spruce"), color = teal,   active = TRUE),
  list(query = intersects, params = list("aspen_birch"), color = teal,   active = TRUE),
  list(query = intersects, params = as.list(allp),        color = orange, active = TRUE))

pdf("Figure_S4.pdf", width = 11, height = 8)
upset(M, sets = allp, nintersects = 6, order.by = "freq",
      mainbar.y.label = "Orthogroup overlap of co-expressologs",
      sets.x.label = "Co-expressolog orthogroups", text.scale = 1.4,
      main.bar.color = "grey35", sets.bar.color = "grey50", queries = q)
dev.off()
cat("Figure_S4.pdf written;", length(allp), "pairwise comparisons\n")
