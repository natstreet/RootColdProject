#!/usr/bin/env Rscript
# external_validation_test.R — enrichment of the conserved-core clusters in published Arabidopsis cold gene sets
# (D062 item 3, D063; design R040 as amended in A035; run on the 63 sample-stable and the 70 sample-sensitive orthogroups
# separately as well as all 133, A036).
#   members     primary: the Col-0 gene of each orthogroup's representative clique (one per orthogroup; kept_lost_orthogroups.tsv);
#               sensitivity: every Col-0 gene in the orthogroup's complete cliques
#   universe    Col-0 genes tested in ComPlEx (in any comparison table with col0)
#   test        hypergeometric, enrichment (upper tail); per stratum and member definition, BH across clusters (1-4 and all
#               clusters pooled) x sets; "match" if q < 0.05 and fold > 1; "not assessable" if expected overlap < 2
#   control     the study's Col-0 cold DEGs (DEGs_col0) must be enriched for calixto_all_de; otherwise stop (mapping broken)
# Inputs (working tree): <sets tsv> (external_gene_sets.R), <kept_lost_orthogroups.tsv>, data/DEGs/DEGs_col0.RData,
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS, data/ComPlEx/PairwiseComparisons/comparison-col0-*-table.RData
# Output (<out dir>): external_validation.tsv, external_validation_control.tsv
# Usage (working tree): Rscript external_validation_test.R <sets tsv> <kept_lost_orthogroups.tsv> <out dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); SETS <- a[1]; KL <- a[2]; OUT <- a[3]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
sets <- fread(SETS); kl <- fread(KL); w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))
uni <- unique(unlist(lapply(list.files("data/ComPlEx/PairwiseComparisons", "^comparison-col0-.*-table\\.RData$", full.names = TRUE),
                            function(f) as.data.table(ld(f))$Species1)))
N <- length(uni)
hyp <- function(mem, set) { mem <- intersect(mem, uni); set <- intersect(set, uni); k <- length(intersect(mem, set)); n <- length(mem); K <- length(set)
  e <- n * K / N; data.table(members = n, set_in_universe = K, overlap = k, expected = e, fold = if (e > 0) k / e else NA_real_,
                             p = phyper(k - 1, K, N - K, n, lower.tail = FALSE)) }
degs <- ld("data/DEGs/DEGs_col0.RData")
ctl <- hyp(degs, sets[set == "calixto_all_de", gene])[, test := "Col-0 DEGs (this study) in calixto_all_de"]
fwrite(ctl, file.path(OUT, "external_validation_control.tsv"), sep = "\t"); print(ctl)
if (!(ctl$p < 0.05 && ctl$fold > 1)) stop("positive control fails: Col-0 DEGs not enriched for the Calixto cold-responsive set")
mem_all <- function(ogs) unique(w[OrthoGroup %in% ogs, col0])
strata <- list(all = kl, stable = kl[core_status == "core"], sensitive = kl[core_status == "sample-sensitive"])
res <- rbindlist(lapply(names(strata), function(st) { s <- strata[[st]]
  rbindlist(lapply(c("representative", "all_clique_members"), function(md) {
    cls <- c(as.list(setNames(1:4, paste0("cluster_", 1:4))), list(all_clusters = 1:4))
    r <- rbindlist(lapply(names(cls), function(cn) { o <- s[cluster %in% cls[[cn]]]
      mem <- if (md == "representative") unique(o$col0_gene) else mem_all(o$OrthoGroup)
      rbindlist(lapply(unique(sets$set), function(se) hyp(mem, sets[set == se, gene])[, `:=`(cluster = cn, set = se)])) }))
    r[, q := p.adjust(p, "BH")][, call := fifelse(expected < 2, "not assessable", fifelse(q < 0.05 & fold > 1, "match", "no match"))]
    r[, `:=`(stratum = st, members_definition = md)] })) }))
setcolorder(res, c("stratum", "members_definition", "cluster", "set")); res[, universe := N]
fwrite(res, file.path(OUT, "external_validation.tsv"), sep = "\t")
print(res[members_definition == "representative", .(stratum, cluster, set, members, overlap, expected = round(expected, 2), fold = round(fold, 2), q = signif(q, 2), call)], nrows = 200)
