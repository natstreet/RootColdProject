#!/usr/bin/env Rscript
# make_TableS2.R — Table S2: GO (biological-process) enrichment of the co-expressolog clade
# sets (the Venn regions of Fig. 4A: ABC / AB / BC / A / B), per species, restricted to that
# species' DEGs.
#
# Enrichment: topGO, **classic** Fisher's exact test, BH adjustment, padj < 0.05.
# Co-expressolog threshold: bidirectional adjusted P < THRESHOLD (default 0.1, the value used in
# the paper). Pass a different value as the first argument.
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/ComPlEx/co_expressologs.RData
#   data/annotation/{Athal_go_ids.tsv, Potra02_go_ids.tsv, Picab02_go_ids.csv}
#   data/annotation/{bg_col0_ATC.rda, bg_aspen_PT.rda, bg_spruce_PA.rda}
#   data/DEGs/DEGs_{col0,ost0,aspen,spruce}.RData
# Output:
#   tables/TableS2_cladeGO_classic_p<thr>.tsv
#
# Usage: Rscript scripts/GO/make_TableS2.R [threshold]     # e.g. 0.1 (default) or 0.05

suppressMessages({ library(topGO); library(dplyr) })

thr <- suppressWarnings(as.numeric(commandArgs(trailingOnly = TRUE)[1])); if (is.na(thr)) thr <- 0.1

load(Sys.glob("data/*/co_expressologs.RData")[1])
mk <- function(a, b) apply(cbind(a, b), 1, function(v) paste(sort(v), collapse = "|"))
co_expressologs$pairkey <- mk(co_expressologs$Species1, co_expressologs$Species2)
ogset <- function(k, p) unique(co_expressologs$OrthologGroup[co_expressologs$pairkey == k & co_expressologs$MaxpVal < p])
regions <- function(p) { A <- ogset("col0|ost0", p); B <- ogset("aspen|birch", p); C <- ogset("pine|spruce", p)
  list(ABC = Reduce(intersect, list(A, B, C)), AB = setdiff(intersect(A, B), C), BC = setdiff(intersect(B, C), A),
       A = setdiff(A, union(B, C)), B = setdiff(B, union(A, C))) }
sp_pair <- list(Col0 = "col0|ost0", Ost0 = "col0|ost0", Aspen = "aspen|birch", Spruce = "pine|spruce")
sp_side <- list(Col0 = "col0", Ost0 = "ost0", Aspen = "aspen", Spruce = "spruce")
coexpr_genes <- function(species, ogs, p) { pk <- sp_pair[[species]]; side <- sp_side[[species]]
  d <- co_expressologs[co_expressologs$pairkey == pk & co_expressologs$MaxpVal < p & co_expressologs$OrthologGroup %in% ogs, ]
  unique(c(d$GeneSpecies1[d$Species1 == side], d$GeneSpecies2[d$Species2 == side])) }

ld <- function(f, o = NULL) { e <- new.env(); load(f, envir = e); if (is.null(o)) get(ls(e)[1], envir = e) else get(o, e) }
athal <- readMappings("data/annotation/Athal_go_ids.tsv")[-1]
potra <- readMappings("data/annotation/Potra02_go_ids.tsv")[-1]
picab_df <- read.csv("data/annotation/Picab02_go_ids.csv"); picab <- unstack(picab_df[, c("go_id", "external_gene_name")])
bg_ATC <- ld("data/annotation/bg_col0_ATC.rda"); bg_PT <- ld("data/annotation/bg_aspen_PT.rda"); bg_PA <- ld("data/annotation/bg_spruce_PA.rda")
ANN <- list(Col0 = athal[names(athal) %in% bg_ATC], Ost0 = athal[names(athal) %in% bg_ATC],
            Aspen = potra[names(potra) %in% bg_PT], Spruce = picab[names(picab) %in% bg_PA])
DEG <- list(Col0 = ld("data/DEGs/DEGs_col0.RData"), Ost0 = ld("data/DEGs/DEGs_ost0.RData"),
            Aspen = ld("data/DEGs/DEGs_aspen.RData"), Spruce = ld("data/DEGs/DEGs_spruce.RData"))

run_go <- function(species, gene_set) { g2 <- ANN[[species]]
  gl <- factor(as.integer(names(g2) %in% gene_set)); names(gl) <- names(g2)
  if (sum(gl == 1) < 3) return(data.frame())
  GOdata <- new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = g2)
  res <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  GenTable(GOdata, classicFisher = res, orderBy = "classicFisher", topNodes = length(usedGO(GOdata))) %>%
    mutate_at("classicFisher", as.numeric) %>% mutate(padj = round(p.adjust(classicFisher, "BH"), 4)) %>%
    arrange(padj) %>% filter(padj < 0.05) }

cells <- list(c("Col0","ABC"), c("Aspen","ABC"), c("Col0","AB"), c("Aspen","AB"),
              c("Aspen","BC"), c("Spruce","BC"), c("Col0","A"), c("Aspen","B"))
reg <- regions(thr); dir.create("tables", showWarnings = FALSE)
all <- list()
cat(sprintf("Table S2 clade GO (classic Fisher), co-expressolog P < %.2f\n", thr))
for (cl in cells) { sp <- cl[1]; rg <- cl[2]
  gs <- intersect(coexpr_genes(sp, reg[[rg]], thr), DEG[[sp]])
  r <- run_go(sp, gs)
  if (nrow(r)) { r$species <- sp; r$set <- rg; all[[paste(sp, rg)]] <- r }
  cat(sprintf("  %-11s study=%-4d sig=%-4d top: %s\n", paste(sp, rg), length(gs), nrow(r), if (nrow(r)) r$Term[1] else "-")) }
out <- do.call(rbind, all)
f <- sprintf("tables/TableS2_cladeGO_classic_p%02d.tsv", round(thr * 100))
write.table(out[, c("species","set","GO.ID","Term","Annotated","Significant","classicFisher","padj")],
            f, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", f, "(", nrow(out), "rows )\n")
