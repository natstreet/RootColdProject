#!/usr/bin/env Rscript
# make_Fig5.R — Figure 5: expression of the conserved-core cliques across species, grouped into
# the four expression clusters. Rows are the 116 complete cliques (one orthogroup each), ordered
# by cluster; columns are the six species x five cold-time-course conditions; values are the
# clique gene's z-scored expression within each species. Row gaps separate the four clusters,
# column gaps separate the species.
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs:
#   data/ComPlEx/cliques/clique_HMcluster_list.rda        clique -> cluster (ID, Cluster)
#   data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS  orthogroup -> per-species gene
#   data/expression/<species>_root_cold_stress_expression.txt
# Output: Figure5.pdf, Figure5.png
# Usage:  Rscript scripts/figures/make_Fig5.R

suppressPackageStartupMessages({ library(pheatmap) })
ld <- function(f){ e<-new.env(); load(f,envir=e); get(ls(e)[1],envir=e) }
species <- c(col0="Col-0", ost0="Ost-0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")
pre     <- c(col0="Col", ost0="Ost", aspen="Pt", birch="Bp", spruce="Pa", pine="Ps")
conds   <- c("01","02","03","04","05")

cc <- ld("data/ComPlEx/cliques/clique_HMcluster_list.rda")     # ID, Cluster
# renumber clusters largest-first to match the manuscript (Cluster 1 = 50, 2 = 30, 3 = 20, 4 = 16)
sz <- sort(table(cc$Cluster), decreasing = TRUE)
cc$Cluster <- as.integer(setNames(seq_along(sz), names(sz))[as.character(cc$Cluster)])
cc <- cc[order(cc$Cluster), ]; cc$OG <- sub("_.*$", "", cc$ID)
cl <- readRDS("data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS")
# orthogroup -> one gene per species
gene_of <- function(og, sp){
  r <- cl[cl$OrthoGroup==og, ]
  g <- unique(c(r$GeneSpecies1[r$Species1==sp], r$GeneSpecies2[r$Species2==sp]))
  if (length(g)) g[1] else NA
}
EXP <- lapply(names(species), function(sp)
  read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp), check.names=FALSE, row.names=1))
names(EXP) <- names(species)
# per gene per species: z-scored 5-condition mean profile
prof <- function(sp, g){
  if (is.na(g) || !(g %in% rownames(EXP[[sp]]))) return(rep(NA, 5))
  v <- as.numeric(EXP[[sp]][g, ]); cc2 <- sub(paste0("^",pre[sp],"\\d+\\."), "", colnames(EXP[[sp]]))
  m <- sapply(conds, function(k) mean(v[cc2==k])); as.numeric(scale(m))
}
M <- t(sapply(seq_len(nrow(cc)), function(i)
  unlist(lapply(names(species), function(sp) prof(sp, gene_of(cc$OG[i], sp))))))
rownames(M) <- cc$ID
colnames(M) <- paste(rep(unname(species), each=5), rep(c("Ctl","6h","24h","3d","10d"), 6), sep=".")
M[M >  2.5] <-  2.5; M[M < -2.5] <- -2.5
M <- M[stats::complete.cases(M), ]

ann <- data.frame(Cluster = factor(cc$Cluster[match(rownames(M), cc$ID)]), row.names = rownames(M))
gaps_row <- head(cumsum(table(ann$Cluster)), -1)
pal <- colorRampPalette(c("#053061","#2166AC","#F7F7F7","#D6604D","#67001F"))(101)
pheatmap(M, cluster_rows=FALSE, cluster_cols=FALSE, gaps_row=gaps_row, gaps_col=c(5,10,15,20,25),
  color=pal, breaks=seq(-2.5,2.5,length.out=102), show_rownames=FALSE,
  annotation_row=ann, border_color=NA, fontsize_col=7,
  main="Conserved-core clique expression clusters",
  filename="Figure5.pdf", width=9, height=10)
pheatmap(M, cluster_rows=FALSE, cluster_cols=FALSE, gaps_row=gaps_row, gaps_col=c(5,10,15,20,25),
  color=pal, breaks=seq(-2.5,2.5,length.out=102), show_rownames=FALSE,
  annotation_row=ann, border_color=NA, fontsize_col=7,
  main="Conserved-core clique expression clusters",
  filename="Figure5.png", width=9, height=10)
cat("cliques plotted:", nrow(M), "| per cluster:", paste(table(ann$Cluster), collapse=" "), "\n")
