#!/usr/bin/env Rscript
# make_TableS4.R — Table S4: GO (biological-process) enrichment of the conserved-core clique
# clusters, per species (manuscript cluster numbering).
#
# The 116 conserved-core cliques (one gene per species, every pair a co-expressolog) are grouped
# into four clusters by the co-expression heatmap (Fig. 5). For each cluster x species the member
# genes are tested for GO enrichment.
#
# Enrichment: topGO, **classic** Fisher's exact test, BH adjustment, padj < 0.05 (consistent with
# Tables S1/S2). Default nodeSize (1). scores = -log10(p). The Table_S4.xlsx p-value column
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
#   tables/TableS4_cliqueClusterGO_classic.tsv
#   tables/Table_S4.xlsx   (the supplementary table: Legend sheet + "Table S4" sheet; padj to 4 s.f., p to 3)
#
# Usage: Rscript scripts/GO/make_TableS4.R

suppressMessages(library(dplyr)); source("scripts/GO/go_table_utils.R")

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
athal <- read_go_mapping(file.path(A, "Athal_go_ids.tsv"))
potra <- read_go_mapping(file.path(A, "Potra02_go_ids.tsv"))
betpen<- read_go_mapping(file.path(A, "Betpen_go_ids.tsv"))
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
  res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  go_result(GOdata, res) %>% mutate(scores = -log10(classicFisher))
}

## manuscript cluster -> analysis cluster
mm <- read.delim("data/ComPlEx/cliques/manuscript_cluster_map.tsv")               # scripts/ComPlEx/manuscript_cluster_map.R
MAP <- setNames(mm$analysis_cluster, mm$manuscript_cluster)                                # manuscript -> analysis
out <- list()
cat("Table S4 clique-cluster GO (classic Fisher, padj<0.05)\n")
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
write.table(D, "tables/TableS4_cliqueClusterGO_classic.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write_supp_xlsx(D, "tables/Table_S4.xlsx", "Table S4", sig_cols = c("classicFisher", "padj", "scores"),
  legend = sprintf(paste("Table S4 GO enrichment of the %s conserved-core clique expression clusters, per species",
                         "(manuscript cluster numbering; topGO classic Fisher's exact test, biological process;",
                         "padj < 0.05).%s"), c("one","two","three","four","five","six","seven","eight","nine")[length(MAP)],
                   {  none <- setdiff(names(MAP), as.character(unique(D[["Cluster (manuscript)"]])))   # clusters without a row
                      if (length(none)) sprintf(" Cluster%s %s had no enriched term.", if (length(none) > 1) "s" else "",
                                                paste(none, collapse = " and ")) else "" }))
cat("Wrote tables/TableS4_cliqueClusterGO_classic.tsv and tables/Table_S4.xlsx (", nrow(D), "rows )\n")
