#!/usr/bin/env Rscript
# hog_clique_consistency.R — for each orthogroup with a complete clique, how many of its complete cliques have all six
# members in one phylogenetic hierarchical orthogroup (HOG) at a given node of the OrthoFinder species tree, and the
# size of that HOG. Assessment for the orthology rule (DECISIONS.md P7): a clique is evidence of conserved
# co-expression between orthologues only if its members descend from one gene of the node's ancestor. This script
# reports; it changes no analysis object. Whether and how the rule enters the core is NRS's decision.
#
# The HOG table must come from the same OrthoFinder run as data/ComPlEx/Orthogroups_20240613.tsv: the script stops
# unless every orthogroup has the same number of genes per species in both tables (as orthogroup_hog_split.R checks
# for N0). HOG gene names are transcripts; they are mapped to the expression matrices' gene identifiers with the
# Phytozome annotation (Arabidopsis PAC ids -> TAIR loci) and by removing transcript suffixes (other species).
#
# Inputs:  data/ComPlEx/Orthogroups_20240613.tsv, a HOG table (N<k>.tsv[.gz]), data/annotation/
#          Athaliana_167_TAIR10.annotation_info.txt.gz, data/expression/<genotype>_root_cold_stress_expression.txt,
#          data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS
# Outputs: <prefix>_per_orthogroup.tsv  orthogroup, expressed_genes, complete_cliques, hog_consistent_cliques,
#                                       hogs (number of HOGs holding its expressed genes), largest_hog_genes,
#                                       consistent_hogs (HOGs holding a consistent clique) and largest_consistent_hog_genes
#                                       (expressed genes of the largest of those: whether size remains a problem)
#          <prefix>_tf_cliques.tsv      for each transcription-factor clique of the current core (tf_cliques.R, at least
#                                       5 of 6 members predicted TFs): its consistent cliques, how many keep at least
#                                       5 TF members, and the members of the best such clique (lowest MeanCliqueSum)
#          <prefix>_summary.tsv         counts
# Usage (repository root): Rscript scripts/analysis/hog_clique_consistency.R <HOG table> <output prefix>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("usage: hog_clique_consistency.R <HOG table> <output prefix>")
hog_file <- a[1]; prefix <- a[2]
rdz <- function(f) if (grepl("[.]gz$", f)) fread(cmd = paste("gzip -dc", shQuote(f))) else fread(f)
SP <- c(Arabidopsis_thaliana = "arabidopsis", Betula_pendula = "birch", Picea_abies = "spruce",
        Pinus_sylvestris = "pine", Populus_tremula = "aspen")
G <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")

h <- rdz(hog_file)
og <- fread("data/ComPlEx/Orthogroups_20240613.tsv", colClasses = "character")
if (!"Ortholog_Group" %in% names(og)) og[[1]] <- NULL
setnames(og, "Ortholog_Group", "OG")
ngenes <- function(v) ifelse(is.na(v) | v == "", 0L, lengths(strsplit(v, ", ", fixed = TRUE)))
sp_all <- intersect(setdiff(names(og), "OG"), names(h))
ca <- og[, c(list(OG = OG), lapply(.SD, ngenes)), .SDcols = sp_all]
cb <- h[, lapply(.SD, function(v) sum(ngenes(v))), by = OG, .SDcols = sp_all]
m <- merge(ca, cb, by = "OG", suffixes = c(".og", ".hog"))
same <- all(sapply(sp_all, function(s) all(m[[paste0(s, ".og")]] == m[[paste0(s, ".hog")]])))
if (!same || nrow(cb) > nrow(ca)) stop("the HOG table is not from the same OrthoFinder run as Orthogroups_20240613.tsv")

pac <- fread(cmd = "gzip -dc data/annotation/Athaliana_167_TAIR10.annotation_info.txt.gz", select = 1:2,
             col.names = c("pac", "locus"), colClasses = "character")
norm <- function(s, g) {
  if (s == "Arabidopsis_thaliana") return(pac$locus[match(sub("^PAC_", "", g), pac$pac)])
  sub("[.][0-9]+$", "", sub("[.]mRNA[.][0-9]+$", "", sub("[.]m[0-9]+$", "", g)))
}
long <- unique(rbindlist(lapply(names(SP), function(s) {
  g <- strsplit(h[[s]], ", ", fixed = TRUE)
  data.table(HOG = rep(h$HOG, lengths(g)), OG = rep(h$OG, lengths(g)), sp = s, gene = norm(s, unlist(g)))
}))[!is.na(gene) & gene != ""])
if (anyDuplicated(long[, .(sp, gene)])) stop("a gene is in more than one HOG")
expr <- unique(unlist(lapply(G, function(x)
  fread(sprintf("data/expression/%s_root_cold_stress_expression.txt", x), select = 1)[[1]])))
e <- long[gene %in% expr]

w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))
hm <- setNames(e$HOG, e$gene)
H <- as.data.table(lapply(G, function(x) unname(hm[w[[x]]]))); setnames(H, G)
unmapped <- H[, Reduce(`|`, lapply(.SD, is.na))]
if (any(unmapped)) stop(sum(unmapped), " complete cliques have a member with no HOG")
w[, one_hog := H[, col0 == ost0 & col0 == aspen & col0 == birch & col0 == spruce & col0 == pine]]
per <- w[, .(complete_cliques = .N, hog_consistent_cliques = sum(one_hog)), by = .(orthogroup = OrthoGroup)]
hs <- e[, .(genes = .N), by = .(OG, HOG)][, .(expressed_genes = sum(genes), hogs = .N, largest_hog_genes = max(genes)),
                                         by = .(orthogroup = OG)]
per <- merge(per, hs, by = "orthogroup", all.x = TRUE)
hsize <- e[, .(hog_genes = .N), by = HOG]
w[, hog := H$col0]
ch <- w[one_hog == TRUE, .(consistent_hogs = uniqueN(hog), largest_consistent_hog_genes = max(hsize$hog_genes[match(unique(hog), hsize$HOG)])),
        by = .(orthogroup = OrthoGroup)]
per <- merge(per, ch, by = "orthogroup", all.x = TRUE)
per[is.na(consistent_hogs), consistent_hogs := 0L]
setcolorder(per, c("orthogroup", "expressed_genes", "complete_cliques", "hog_consistent_cliques", "hogs", "largest_hog_genes",
                   "consistent_hogs", "largest_consistent_hog_genes"))
setorder(per, -expressed_genes, orthogroup)
fwrite(per, paste0(prefix, "_per_orthogroup.tsv"), sep = "\t")

# transcription-factor cliques: the same predicted-TF sets as tf_cliques.R
source("scripts/analysis/tf_cliques.R")
tfc <- as.data.table(tf_cliques())
strip <- function(x) sub("[.][0-9]+$", "", as.character(x))
TF <- list()
for (s in c("col0", "aspen", "birch", "spruce", "pine")) {
  en <- new.env(); load(sprintf("data/TranscriptionFactors/TF_%s.RData", s), envir = en); x <- get(ls(en)[1], en)
  TF[[s]] <- unique(c(strip(x[[1]]), strip(x[[2]])))
}
TF$ost0 <- TF$col0
w[, n_tf := Reduce(`+`, lapply(G, function(x) strip(w[[x]]) %in% TF[[x]]))]
tfr <- rbindlist(lapply(seq_len(nrow(tfc)), function(k) {
  o <- tfc$OrthoGroup[k]; x <- w[OrthoGroup == o]
  ok <- x[one_hog == TRUE & n_tf >= 5][order(MeanCliqueSum)]
  data.table(tf = tfc$tf[k], orthogroup = o, complete_cliques = nrow(x), hog_consistent_cliques = sum(x$one_hog),
             consistent_with_5_tf = nrow(ok),
             best_members = if (nrow(ok)) paste(unlist(ok[1, ..G]), collapse = ", ") else "")
}))
fwrite(tfr, paste0(prefix, "_tf_cliques.tsv"), sep = "\t")
print(tfr[, .(tf, orthogroup, complete_cliques, hog_consistent_cliques, consistent_with_5_tf)], row.names = FALSE)
summ <- data.table(quantity = c("hog_table", "expressed_genes_with_a_hog", "expressed_genes", "orthogroups_with_complete_clique",
                                "orthogroups_with_hog_consistent_clique", "complete_cliques", "hog_consistent_cliques"),
                   value = c(basename(hog_file), sum(expr %in% long$gene), length(expr), nrow(per),
                             sum(per$hog_consistent_cliques > 0), nrow(w), sum(w$one_hog)))
fwrite(summ, paste0(prefix, "_summary.tsv"), sep = "\t")
print(summ, row.names = FALSE)
