#!/usr/bin/env Rscript
# make_FigS1.R — Figure S1: correlation matrix of the 36 per-species expression clusters of
# cold-responsive transcription factors (six clusters x six genotypes), grouped into the twelve
# superclusters.
#
# Each cluster profile is the mean VST expression of its member genes at each timepoint. The matrix
# is ordered by hierarchical clustering (Ward's method, ward.D2, on 1 - Pearson r) and the tree is
# cut into twelve groups, which are the superclusters. Superclusters are numbered 1-12 in dendrogram
# order, the numbering used by the supercluster gene memberships (data/superclusters/membership/);
# the four discussed in the main text are named from SC_mapping.csv (SC_I-SC_IV).
#
# Inputs (run from a directory containing the unpacked figshare deposit `data/`, into which this
# repo's committed data/ files merge):
#   data/superclusters/merged_super_clusters_means.RData   36 x 5 matrix of cluster-mean profiles
#   data/superclusters/SC_mapping.csv                      main-text name -> supercluster number
# Output: figures/Figure_S1.pdf, figures/Figure_S1.png, and the printed supercluster composition.
#
# Usage: Rscript scripts/figures/make_FigS1.R

suppressMessages(library(corrplot))

ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
means <- as.matrix(ld("data/superclusters/merged_super_clusters_means.RData"))   # 36 x 5
lab <- function(x) { x <- sub("_birch$", "_Birch", x); x <- sub("_poplar$", "_Aspen", x); sub("_pine$", "_Pine", x) }
rownames(means) <- lab(rownames(means))
M <- cor(t(means), method = "pearson")

hc  <- hclust(as.dist(1 - M), method = "ward.D2")
ord <- hc$order
grp <- cutree(hc, k = 12)[ord]
sc  <- match(grp, unique(grp))                       # supercluster number, 1-12 in dendrogram order
map <- read.csv("data/superclusters/SC_mapping.csv", stringsAsFactors = FALSE)
main_name <- setNames(map$Main_text, as.integer(sub("^SC", "", map$Original_supercluster)))

for (k in 1:12) cat(sprintf("Supercluster %-2d %-7s %s\n", k,
  if (!is.na(main_name[as.character(k)])) paste0("(", main_name[as.character(k)], ")") else "",
  paste(rownames(M)[ord][sc == k], collapse = ", ")))

n <- nrow(M)
draw <- function() {
  corrplot(M[ord, ord], method = "square", order = "original", type = "lower", outline = TRUE,
           addgrid.col = "darkgray", tl.pos = "l", tl.col = "black", tl.cex = 0.7, cl.ratio = 0.12,
           cl.pos = "b", col = colorRampPalette(c("darkred", "white", "midnightblue"))(100),
           col.lim = c(-1, 1), mar = c(1, 0, 1, 10))
  par(xpd = NA)
  for (k in 1:12) {                                  # supercluster boxes and labels
    r <- range(which(sc == k))
    rect(r[1] - 0.5, n - r[2] + 0.5, r[2] + 0.5, n - r[1] + 1.5, border = "black", lwd = 3)
    nm <- main_name[as.character(k)]
    text(r[2] + 0.8, n + 1 - mean(r), adj = 0, cex = 0.8, font = 2,
         labels = paste0("Supercluster ", k, if (!is.na(nm)) paste0(" (", nm, ")") else ""))
  }
}
dir.create("figures", showWarnings = FALSE)
pdf("figures/Figure_S1.pdf", width = 11, height = 10); draw(); invisible(dev.off())
png("figures/Figure_S1.png", width = 11, height = 10, units = "in", res = 150); draw(); invisible(dev.off())
