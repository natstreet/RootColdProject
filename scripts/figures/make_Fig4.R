#!/usr/bin/env Rscript
# make_Fig4.R — Figure 4: heatmap of the conserved-core cliques across species.
# For each clique, the per-gene z-scored MEDIAN expression of its member in each species over the
# cold time course (Control, 6h, 24h, 3d, 10d). Rows are clustered (Euclidean distance, Ward.D2)
# and cut into four clusters, as in the original analysis (E. van Zalen). The figure shows one block per
# manuscript cluster, in the order the dendrogram places them, joined by one dendrogram; colours are capped at
# |z| = 1.5.
#
# Each block is captioned with GO terms listed by ID per manuscript cluster in a caption table (default
# scripts/figures/Fig4_captions.tsv; fig4_caption_candidates.R proposes terms). Term names are read from the
# Table S3 output, and the script stops unless each term is significant in that cluster, so run
# scripts/GO/make_TableS3.R first.
#
# Inputs:
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS  the core cliques, one gene per species
#   data/ComPlEx/cliques/clique_HMcluster_list.rda              clique -> analysis cluster (clique_clusters.R)
#   scripts/ComPlEx/clique_matrix.R                            the z-scored matrix and the row clustering
#   data/expression/med_expression/<sp>_med.txt.gz             per-timepoint median expression
#   data/ComPlEx/cliques/manuscript_cluster_map.tsv            manuscript -> analysis cluster map
#   tables/TableS3_cliqueClusterGO_classic.tsv                 GO terms per manuscript cluster
#   scripts/figures/Fig4_captions.tsv (or the first argument)   GO IDs captioning each cluster (cluster, GO.ID)
# Output: Figure4.pdf, Figure4.png
# Usage:  Rscript scripts/figures/make_Fig4.R [caption table]

suppressPackageStartupMessages({ library(dplyr); library(readr); library(pheatmap); library(cowplot)
  library(grid); library(ggplot2) })

source("scripts/ComPlEx/clique_matrix.R")                                   # the matrix, shared with clique_clusters.R
sp <- clique_sp; labs <- clique_labs; conds <- clique_conds
ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

## ---- z-scored median expression, one row per clique ----
M  <- clique_matrix()
cat("Cliques in heatmap:", nrow(M), "\n")

## ---- row clustering; check it against the deposited cluster assignment ----
tr  <- clique_tree(M)
ord <- tr$order
cut <- cutree(tr, k = 4)
hm  <- ld(Sys.glob("data/*/cliques/clique_HMcluster_list.rda")[1])            # ID, Cluster (analysis)
ana <- hm$Cluster[match(rownames(M), hm$ID)]
tab <- table(cut, ana)                                                        # must be a one-to-one relabelling
if (anyNA(ana) || any(rowSums(tab > 0) != 1) || any(colSums(tab > 0) != 1))
  stop("row clustering differs from clique_HMcluster_list.rda")
mm <- read.delim("data/ComPlEx/cliques/manuscript_cluster_map.tsv")               # scripts/ComPlEx/manuscript_cluster_map.R
MAP <- setNames(mm$analysis_cluster, mm$manuscript_cluster)                                # manuscript -> analysis
man <- as.integer(names(MAP))[match(ana, MAP)]                               # manuscript cluster per row
mo  <- man[ord]                                                               # manuscript cluster, dendrogram order
runs <- rle(mo)
if (anyDuplicated(runs$values)) stop("a cluster is not contiguous in dendrogram order")
blk_cl <- runs$values                                                         # one block per cluster, top to bottom
bo  <- rep(seq_along(runs$values), runs$lengths)
nb  <- runs$lengths; NB <- length(nb)
cat("Blocks (cluster: size):", paste(blk_cl, nb, sep = ": ", collapse = ", "), "\n")

## ---- block captions: GO IDs -> names from Table S3, each required to be significant there ----
s3  <- read.delim("tables/TableS3_cliqueClusterGO_classic.tsv", check.names = FALSE)
cargs <- commandArgs(trailingOnly = TRUE)
cap <- read.delim(if (length(cargs)) cargs[1] else "scripts/figures/Fig4_captions.tsv", stringsAsFactors = FALSE)
go_label <- sapply(blk_cl, function(mc) {
  ids <- cap$GO.ID[cap$cluster == mc]
  hit <- s3[s3[["Cluster (manuscript)"]] == mc & s3$GO.ID %in% ids, ]
  miss <- setdiff(ids, hit$GO.ID)
  if (length(miss)) stop("caption table: not significant in Table S3 for cluster ", mc, ": ", paste(miss, collapse = ", "))
  nm <- hit$Term[match(ids, hit$GO.ID)]
  paste(paste0(toupper(substr(nm, 1, 1)), substring(nm, 2)), collapse = "\n") })
tag <- paste("Cluster", blk_cl)

## ---- drawing ----
pal  <- colorRampPalette(c("#313695", "#f8f8f8", "#A50026"))(100)
brks <- seq(-1.5, 1.5, length.out = 101)
Mo   <- pmin(pmax(M[ord, ], -1.5), 1.5)                                      # rows in dendrogram order
spl  <- split(seq_len(nrow(Mo)), bo)
mkblock <- function(rows) pheatmap(Mo[rows, , drop = FALSE], cluster_rows = FALSE, cluster_cols = FALSE,
  gaps_col = c(5, 10, 15, 20, 25), color = pal, breaks = brks, legend = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, border_color = NA, silent = TRUE)$gtable
gb <- lapply(spl, mkblock)
leg <- pheatmap(Mo, cluster_rows = FALSE, cluster_cols = FALSE, color = pal, breaks = brks,
  legend = TRUE, show_rownames = FALSE, show_colnames = FALSE, border_color = NA, silent = TRUE)$gtable
leg <- leg$grobs[[which(leg$layout$name == "legend")]]

CW <- 9.5; CH <- 10.6
BX <- 0.085; BW <- 0.805; Wreg <- BW * CW                                    # matrix-block region on canvas
L  <- gb[[1]]$layout; mi <- which(L$name == "matrix")
absin <- sapply(seq_along(gb[[1]]$widths), function(k) convertWidth(gb[[1]]$widths[[k]], "in", valueOnly = TRUE))
absin[L$l[mi]:L$r[mi]] <- 0
mx0 <- sum(absin[seq_len(L$l[mi] - 1)]) / Wreg; mx1 <- mx0 + (Wreg - sum(absin)) / Wreg
mL <- BX + BW * mx0; mR <- BX + BW * mx1
xm <- function(fr) BX + BW * (mx0 + (mx1 - mx0) * fr)
gf <- 0.24; tot_u <- 30 + 5 * gf
tp_cen <- sapply(1:30, function(c){ g <- sum(c > c(5, 10, 15, 20, 25)); (c - 0.5 + g * gf) / tot_u })
sp_cen <- sapply(1:6,  function(k) (((k - 1) * 5 + (k - 1) * gf) + 2.5) / tot_u)

hH <- 0.045; hG <- 0.058; hT <- 0.052; hS <- 0.030; topm <- 0.008; botm <- 0.006
avail <- 1 - topm - hH - (NB - 1) * hG - hT - hS - botm
bh <- avail * nb / sum(nb)
y_top <- (1 - topm - hH) - c(0, cumsum(bh[-NB] + hG))
y_bot <- y_top - bh
y_cap <- c((1 - topm + y_top[1]) / 2, y_bot[-NB] - hG / 2)
y_tp <- y_bot[NB] - hT / 2; y_sp <- y_bot[NB] - hT - hS / 2

## one dendrogram across the blocks: leaf position -> canvas y, gaps built in
yj <- unlist(lapply(seq_len(NB), function(b) y_top[b] - (seq_len(nb[b]) - 0.5) / nb[b] * bh[b]))
ymap <- approxfun(seq_along(yj), yj, rule = 2)
lp <- match(seq_len(nrow(M)), tr$order); nodex <- numeric(nrow(tr$merge)); H <- tr$height
gx <- function(k) if (k < 0) lp[-k] else nodex[k]; gy <- function(k) if (k < 0) 0 else H[k]
seg <- NULL
for (i in seq_len(nrow(tr$merge))) {
  a <- tr$merge[i, 1]; b <- tr$merge[i, 2]; xa <- gx(a); xb <- gx(b); h <- H[i]; nodex[i] <- (xa + xb) / 2
  seg <- rbind(seg, cbind(x = c(xa, xb, xa), y = c(gy(a), gy(b), h), xe = c(xa, xb, xb), ye = c(h, h, h)))
}
xdR <- mL - 0.006; xdL <- 0.006; xd <- function(h) xdR - (h / max(H)) * (xdR - xdL)
segdf <- data.frame(x = xd(seg[, "y"]), y = ymap(seg[, "x"]), xend = xd(seg[, "ye"]), yend = ymap(seg[, "xe"]))

p <- ggdraw() + theme(plot.background = element_rect(fill = "white", colour = NA))
for (b in seq_len(NB)) p <- p + draw_grob(gb[[b]], x = BX, y = y_bot[b], width = BW, height = bh[b]) +
  draw_label(go_label[b], x = mL + 0.004, y = y_cap[b], hjust = 0, vjust = 0.5, size = 11, fontface = "bold") +
  draw_label(tag[b], x = mR - 0.004, y = y_cap[b], hjust = 1, vjust = 0.5, size = 11.5, fontface = "bold", colour = "#8C8C8C")
p <- p + draw_grob(leg, x = 0.905, y = 0.34, width = 0.09, height = 0.34) +
  geom_segment(data = segdf, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.22, colour = "grey35")
for (c in 1:30) p <- p + draw_label(rep(conds, 6)[c], x = xm(tp_cen[c]), y = y_tp, angle = 90, hjust = 1, vjust = 0.5, size = 8)
for (k in 1:6)  p <- p + draw_label(labs[[k]], x = xm(sp_cen[k]), y = y_sp, hjust = 0.5, vjust = 0.5, size = 10, fontface = "bold")

pdf("Figure4.pdf", width = CW, height = CH, encoding = "ISOLatin1"); print(p); invisible(dev.off())
ggsave("Figure4.png", p, width = CW, height = CH, dpi = 300, bg = "white")
