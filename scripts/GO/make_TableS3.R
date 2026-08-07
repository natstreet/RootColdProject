#!/usr/bin/env Rscript
# make_TableS3.R — Table S3: GO (biological-process) enrichment of the conserved-core clique
# clusters, per species (manuscript cluster numbering).
#
# The 116 conserved-core cliques (one gene per species, every pair a co-expressolog) are grouped
# into four clusters by the co-expression heatmap (Fig. 5). For each cluster x species the member
# genes are tested for GO enrichment.
#
# Enrichment: topGO, **classic** Fisher's exact test, BH adjustment, padj < 0.05 (consistent with
# Tables S1/S2). Default nodeSize (1). scores = -log10(p). The Table_S3.xlsx p-value column
# (classicFisher) holds these classic-Fisher values; regenerated numbers match the deposited table
# to within GO.db-version rounding (e.g. GO:0008283 4.2e-10, 59 terms).
#
# Cluster numbering: the heatmap/analysis clusters are renumbered for the manuscript
#   manuscript 1 = analysis 2 (cell division) ; 2 = analysis 4 (metabolic) ;
#   3 = analysis 3 (signalling)               ; 4 = analysis 1 (16 cell-wall OGs, no enrichment).
# Cluster 4 yields no enriched terms and so does not appear in the table.
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/ComPlEx/cliques/clique_HMcluster_list.rda              (clique_cluster_order: cliqueID -> analysis cluster)
#   data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS   (per-species gene of each clique)
#   data/annotation/{Athal,Potra02,Betpen,Picab02,Pinsy01}_go_ids.{tsv,csv}
#   data/annotation/bg_{col0_ATC,ost0_ATO,aspen_PT,birch_BP,spruce_PA,pine_PS}.rda
# Output:
#   tables/TableS3_cliqueClusterGO_classic.tsv
#
# Usage: Rscript scripts/GO/make_TableS3.R

suppressMessages({ library(topGO); library(dplyr) })

ld <- function(f, o = NULL) { e <- new.env(); load(f, envir = e); if (is.null(o)) get(ls(e)[1], envir = e) else get(o, e) }

## --- clique cluster membership -> per (analysis cluster, species) gene sets ---
cc  <- ld(Sys.glob("data/*/cliques/clique_HMcluster_list.rda")[1])          # ID, Cluster (analysis)
cg  <- readRDS(Sys.glob("data/*/cliques/clique_genes_filterable_COMPLETE.RDS")[1])
# long pairwise table -> one gene per (cliqueID, species)
glong <- unique(rbind(
  data.frame(cliqueID = cg$cliqueID, sp = cg$Species1, gene = cg$GeneSpecies1, stringsAsFactors = FALSE),
  data.frame(cliqueID = cg$cliqueID, sp = cg$Species2, gene = cg$GeneSpecies2, stringsAsFactors = FALSE)))
glong$Cluster <- cc$Cluster[match(glong$cliqueID, cc$ID)]
glong <- glong[!is.na(glong$Cluster), ]
gene_set <- function(anaclust, sp) unique(glong$gene[glong$Cluster == anaclust & glong$sp == sp])

## --- per-species GO annotation restricted to the enrichment background ---
A <- "data/annotation"
athal <- readMappings(file.path(A, "Athal_go_ids.tsv"))[-1]
potra <- readMappings(file.path(A, "Potra02_go_ids.tsv"))[-1]
betpen<- readMappings(file.path(A, "Betpen_go_ids.tsv"))[-1]
picab <- unstack(read.csv(file.path(A, "Picab02_go_ids.csv"))[, c("go_id","external_gene_name")])
pinsy <- unstack(read.csv(file.path(A, "Pinsy01_go_ids.csv"))[, c("go_id","external_gene_name")])
bg <- function(o) ld(file.path(A, sprintf("bg_%s.rda", o)))
ANN <- list(
  col0   = athal [names(athal)  %in% bg("col0_ATC")],
  ost0   = athal [names(athal)  %in% bg("ost0_ATO")],
  aspen  = potra [names(potra)  %in% bg("aspen_PT")],
  birch  = betpen[names(betpen) %in% bg("birch_BP")],
  spruce = picab [names(picab)  %in% bg("spruce_PA")],
  pine   = pinsy [names(pinsy)  %in% bg("pine_PS")])
LAB <- c(col0="Col0", ost0="Ost0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")

run_go <- function(sp, genes) {
  g2 <- ANN[[sp]]
  gl <- factor(as.integer(names(g2) %in% genes)); names(gl) <- names(g2)
  if (sum(gl == 1) < 1) return(data.frame())
  GOdata <- new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = g2)
  allGO <- usedGO(GOdata)
  res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  GenTable(GOdata, classicFisher = res, orderBy = "classicFisher", topNodes = length(allGO)) %>%
    mutate_at("classicFisher", as.numeric) %>%
    mutate(padj = round(p.adjust(classicFisher, "BH"), 4)) %>%
    arrange(padj) %>% filter(padj < 0.05) %>%
    mutate(scores = -log10(classicFisher))
}

## manuscript cluster -> analysis cluster
MAP <- c(`1` = 2, `2` = 4, `3` = 3, `4` = 1)
out <- list()
cat("Table S3 clique-cluster GO (classic Fisher, padj<0.05)\n")
for (mc in names(MAP)) for (sp in names(LAB)) {
  gs <- gene_set(MAP[[mc]], sp)
  r <- run_go(sp, gs)
  if (nrow(r)) {
    r <- transform(r, `Cluster (manuscript)` = as.integer(mc), Species = LAB[[sp]],
                   analysis_cluster = MAP[[mc]], check.names = FALSE)
    out[[paste(mc, sp)]] <- r
  }
  cat(sprintf("  C%s(ana%d) %-7s study=%-3d sig=%-3d\n", mc, MAP[[mc]], sp, length(gs), nrow(r)))
}
D <- do.call(rbind, out)
D <- D[, c("Cluster (manuscript)","Species","GO.ID","Term","Annotated","Significant","Expected",
           "classicFisher","padj","scores","analysis_cluster")]
dir.create("tables", showWarnings = FALSE)
write.table(D, "tables/TableS3_cliqueClusterGO_classic.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote tables/TableS3_cliqueClusterGO_classic.tsv (", nrow(D), "rows )\n")
