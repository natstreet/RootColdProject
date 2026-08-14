#!/usr/bin/env Rscript
# make_FigS1.R — Figure S1: correlation matrix of the 36 per-species expression clusters
# (six clusters x six genotypes).
#
# The correlation matrix is fully reproducible from the deposited cluster-mean profiles (below).
# The twelve cross-species superclusters shown as boxes on the published figure are a CURATED
# grouping of these 36 clusters: it is NOT recoverable automatically (neither an hclust cutree at
# k=12 nor maximum gene-list overlap reproduces it, because the supercluster gene sets overlap
# heavily). The grouping and its labels are therefore a manual annotation layer; the assignment is
# documented by the deposited supercluster memberships and `SC_mapping.csv`
# (main-text SC_I-SC_IV = superclusters 6, 2, 3, 1). If a curated
# `superclusters/cluster_to_supercluster.csv` (columns: cluster, supercluster) is present, the boxes
# are drawn from it; otherwise the bare matrix is produced.
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/superclusters/merged_super_clusters_means.RData   (36 x 5 matrix: cluster x timepoint)
#   data/superclusters/cluster_to_supercluster.csv         (optional; curated grouping for the boxes)
# Output:
#   figures/Figure_S1.pdf
#
# Usage: Rscript scripts/figures/make_FigS1.R

suppressMessages(library(corrplot))

ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
means <- ld(Sys.glob("data/superclusters/merged_super_clusters_means.RData")[1])   # 36 x 5
M <- cor(t(means), method = "pearson")

grp_file <- "data/superclusters/cluster_to_supercluster.csv"
have_grp <- file.exists(grp_file)
if (have_grp) {
  g <- read.csv(grp_file, stringsAsFactors = FALSE)
  g <- g[match(rownames(M), g$cluster), ]
  ord <- order(g$supercluster); M <- M[ord, ord]
  rect_n <- as.integer(table(factor(g$supercluster[ord], levels = unique(g$supercluster[ord]))))
}

dir.create("figures", showWarnings = FALSE)
pdf("figures/Figure_S1.pdf", width = 11, height = 9)
corrplot(M, method = "square", type = "upper", order = if (have_grp) "original" else "hclust",
         tl.col = "black", tl.cex = 0.7, tl.srt = 45,
         col = colorRampPalette(c("#7F0000", "white", "#00007F"))(200), cl.pos = "r")
if (have_grp) corrRect(rect_n, col = "black", lwd = 2)
dev.off()
cat("Wrote figures/Figure_S1.pdf",
    if (have_grp) "(with curated supercluster boxes)" else "(matrix only; supply cluster_to_supercluster.csv for boxes)", "\n")
