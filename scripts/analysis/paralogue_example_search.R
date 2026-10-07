#!/usr/bin/env Rscript
# paralogue_example_search.R — the bounded search for an illustration of a paralogue-carried core orthogroup (D045),
# with criteria fixed before the search (ledger R035) and not tuned afterwards.
#
# Candidates: the paralogue-carried core orthogroups with <= 20 expressed genes (results/paralogue_carried_small.tsv).
# Tree: the IQ-TREE alignment tree (UFBoot support on each branch), rooted on the bipartition at the root of OrthoFinder's
# resolved tree for the same orthogroup. If that bipartition is not a branch of the IQ-TREE tree, the orthogroup cannot be
# assessed as written and is reported so (no other rooting is substituted).
# The node: where the representative clique's members meet (their most recent common ancestor). It must separate them,
# the angiosperm members (Col-0, Ost-0, aspen, birch) wholly in one child clade and the conifer members (spruce, pine)
# wholly in another. A node's support is the UFBoot of the branch above it; for the root, the branch rooted on.
# Qualifies if all of:
#   (a) UFBoot >= 95 on that node;
#   (b) species overlap among the seed-plant genes on both sides: at least two gymnosperm species, or at least two
#       angiosperm species, present on both sides (so one stray gene cannot fake a duplication);
#   (c) the angiosperm members and the conifer members each lie in a clade with UFBoot >= 95 on their own side (the best
#       support among the clades from the group's own common ancestor up to the child of the node is reported).
# More than one qualifies: the most DE genotypes, then the fewest expressed genes. None: reported as none.
#
# Inputs (working directory): results/paralogue_carried_small.tsv, results/orthology_class.tsv,
#   tables/core_response_patterns.tsv, data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS,
#   data/annotation/Athaliana_167_TAIR10.annotation_info.txt.gz, data/orthofinder/N1_msa.tsv
#   <msa tree dir>/<OG>_tree.txt (IQ-TREE), <resolved tree dir>/<OG>_tree.txt (OrthoFinder Resolved_Gene_Trees)
# Outputs: <out>/paralogue_example_candidates.tsv (one row per candidate, each criterion with its values),
#          <out>/paralogue_example_chosen.tsv (the chosen orthogroup's members and N1 HOGs; empty if none),
#          <out>/<OG>_rooted_tree.pdf for the chosen orthogroup
# Usage: Rscript scripts/analysis/paralogue_example_search.R <msa tree dir> <resolved tree dir> <out dir>
suppressPackageStartupMessages({ library(ape); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); MSA <- a[1]; RES <- a[2]; OUT <- a[3]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
UF <- 95; MIN_SP <- 2
ANGIO <- c("Amborella_trichopoda","Arabidopsis_thaliana","Beta_vulgaris","Betula_pendula","Camellia_sinensis",
           "Cinnamomum_kanehirae","Cucumis_sativus","Eucalyptus_grandis","Malania_oleifera","Malus_domestica",
           "Medicago_truncatula","Mimulus_guttatus","Nicotiana_tabacum","Populus_tremula","Prunus_avium","Salix_purpurea",
           "Vitis_vinifera","Zostera_marina")
GYMNO <- c("Ginkgo_biloba","Picea_abies","Pinus_sylvestris","Pinus_tabuliformis","Sequoiadendron_giganteum",
           "Taxus_chinensis","Torreya_grandis","Welwitschia_mirabilis")
STUDY <- c(col0 = "Arabidopsis_thaliana", ost0 = "Arabidopsis_thaliana", aspen = "Populus_tremula",
           birch = "Betula_pendula", spruce = "Picea_abies", pine = "Pinus_sylvestris")
SPP <- c(ANGIO, GYMNO, "Physcomitrella_patens")
sp_of <- function(t) vapply(t, function(x) { s <- SPP[startsWith(x, paste0(SPP, "_"))]; if (length(s)) s[which.max(nchar(s))] else NA_character_ }, "")
pac <- fread(cmd = "gzip -dc data/annotation/Athaliana_167_TAIR10.annotation_info.txt.gz", select = 1:2,
             col.names = c("pac", "locus"), colClasses = "character")
norm <- function(s, g) if (s == "Arabidopsis_thaliana") pac$locus[match(sub("^PAC_", "", g), pac$pac)] else
  sub("[.][0-9]+$", "", sub("[.]mRNA[.][0-9]+$", "", sub("[.]m[0-9]+$", "", g)))
strip <- function(x) sub("[.][0-9]+$", "", x)

cand <- fread("results/paralogue_carried_small.tsv")
core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
resp <- fread("tables/core_response_patterns.tsv")[, .(OrthoGroup, n_genotypes_DE, DE_genotypes)]
n1 <- fread("data/orthofinder/N1_msa.tsv")

clade <- function(tr, k) if (k <= Ntip(tr)) tr$tip.label[k] else extract.clade(tr, k)$tip.label
supp <- function(tr, k) { if (k <= Ntip(tr)) return(NA_real_); suppressWarnings(as.numeric(tr$node.label[k - Ntip(tr)])) }
split_key <- function(tips, all) { s <- sort(tips); o <- sort(setdiff(all, tips)); if (s[1] < o[1]) paste(s, collapse = "|") else paste(o, collapse = "|") }

assess <- function(og) {
  row <- list(OrthoGroup = og)
  it <- read.tree(file.path(MSA, paste0(og, "_tree.txt"))); rt <- read.tree(file.path(RES, paste0(og, "_tree.txt")))
  if (!setequal(it$tip.label, rt$tip.label)) return(c(row, list(status = "tip sets differ between the two trees")))
  # the root bipartition of OrthoFinder's resolved tree, located among the IQ-TREE tree's branches
  side <- clade(rt, rt$edge[rt$edge[, 1] == Ntip(rt) + 1, 2][1])
  ut <- unroot(it); key <- split_key(side, ut$tip.label)
  e <- which(vapply(seq_len(nrow(ut$edge)), function(i) split_key(clade(ut, ut$edge[i, 2]), ut$tip.label) == key, TRUE))
  if (!length(e)) return(c(row, list(status = "OrthoFinder's root bipartition is not a branch of the IQ-TREE tree")))
  root_support <- if (ut$edge[e[1], 2] > Ntip(ut)) supp(ut, ut$edge[e[1], 2]) else NA_real_
  tr <- root(it, outgroup = side, resolve.root = TRUE, edgelabel = TRUE)
  # the representative clique's members as tips
  cl <- core[OrthoGroup == og]; sp <- sp_of(tr$tip.label)
  gene <- mapply(function(t, s) if (is.na(s)) NA_character_ else norm(s, substring(t, nchar(s) + 2)), tr$tip.label, sp)
  tip_of <- function(x) tr$tip.label[which(sp == STUDY[[x]] & gene == strip(cl[[x]]))]
  ang <- unique(unlist(lapply(c("col0", "ost0", "aspen", "birch"), tip_of))); con <- unique(unlist(lapply(c("spruce", "pine"), tip_of)))
  if (length(con) != 2 || length(ang) < 3) return(c(row, list(status = "clique members not all found in the tree")))
  X <- getMRCA(tr, c(ang, con)); kids <- tr$edge[tr$edge[, 1] == X, 2]
  ka <- kids[vapply(kids, function(k) all(ang %in% clade(tr, k)), TRUE)]; kc <- kids[vapply(kids, function(k) all(con %in% clade(tr, k)), TRUE)]
  separated <- length(ka) == 1 && length(kc) == 1 && ka != kc
  geno_tip <- unlist(lapply(names(STUDY), function(x) setNames(tip_of(x), rep(x, length(tip_of(x))))))
  by_side <- paste(vapply(seq_along(kids), function(i) paste(unique(names(geno_tip)[geno_tip %in% clade(tr, kids[i])]), collapse = ","), ""),
                   collapse = " | ")
  sX <- if (X == Ntip(tr) + 1) root_support else supp(tr, X)
  seed <- function(t) { s <- sp_of(t); s[!is.na(s) & s != "Physcomitrella_patens"] }
  if (separated) {
    s1 <- unique(seed(clade(tr, ka))); s2 <- unique(seed(unlist(lapply(setdiff(kids, ka), function(k) clade(tr, k)))))
    both <- intersect(s1, s2); g_both <- sum(both %in% GYMNO); a_both <- sum(both %in% ANGIO)
    best_up <- function(mem, top) { k <- if (length(mem) == 1) match(mem, tr$tip.label) else getMRCA(tr, mem); best <- NA_real_
      repeat { v <- supp(tr, k); if (!is.na(v)) best <- max(best, v, na.rm = TRUE); if (k == top) break; k <- tr$edge[tr$edge[, 2] == k, 1] }
      best }
    ua <- best_up(ang, ka); uc <- best_up(con, kc)
  } else { g_both <- a_both <- NA_integer_; ua <- uc <- NA_real_ }
  pa <- !is.na(sX) && sX >= UF; pb <- separated && (g_both >= MIN_SP || a_both >= MIN_SP)
  pc <- separated && !is.na(ua) && !is.na(uc) && ua >= UF && uc >= UF
  c(row, list(status = "assessed", tips = Ntip(tr), node_is_root = X == Ntip(tr) + 1, members_by_side = by_side,
              members_separated = separated,
              node_ufboot = sX, a_pass = pa, gymnosperm_species_both_sides = g_both, angiosperm_species_both_sides = a_both,
              b_pass = pb, angiosperm_members_best_ufboot = ua, conifer_members_best_ufboot = uc, c_pass = pc,
              qualifies = separated && pa && pb && pc))
}
res <- rbindlist(lapply(cand$OrthoGroup, assess), fill = TRUE)
res <- merge(cand[, .(OrthoGroup, expressed_genes, col0_symbol)], res, by = "OrthoGroup")
res <- merge(res, resp, by = "OrthoGroup", all.x = TRUE)
res[is.na(qualifies), qualifies := FALSE]
setorder(res, -qualifies, -n_genotypes_DE, expressed_genes, OrthoGroup)
fwrite(res, file.path(OUT, "paralogue_example_candidates.tsv"), sep = "\t", na = "NA")
print(res[, .(OrthoGroup, expressed_genes, status, node_is_root, members_by_side, members_separated, node_ufboot, a_pass,
              gymnosperm_species_both_sides, angiosperm_species_both_sides, b_pass, angiosperm_members_best_ufboot,
              conifer_members_best_ufboot, c_pass, qualifies, n_genotypes_DE)], row.names = FALSE)

q <- res[qualifies == TRUE]
if (!nrow(q)) {
  writeLines("no candidate qualifies", file.path(OUT, "paralogue_example_chosen.tsv"))
  cat("no candidate qualifies\n")
} else {
  og <- q$OrthoGroup[1]; cl <- core[OrthoGroup == og]
  hog <- function(s, g) { r <- n1[OG == og]; for (k in seq_len(nrow(r))) { v <- strsplit(r[[s]][k], ", ", fixed = TRUE)[[1]]
    if (length(v) && strip(g) %in% vapply(v, function(t) norm(s, t), "")) return(r$HOG[k]) }; NA_character_ }
  ch <- data.table(OrthoGroup = og, genotype = names(STUDY), gene = unlist(cl[, names(STUDY), with = FALSE]))
  ch[, N1_HOG := mapply(function(x, g) hog(STUDY[[x]], g), genotype, gene)]
  fwrite(ch, file.path(OUT, "paralogue_example_chosen.tsv"), sep = "\t")
  print(ch, row.names = FALSE)
  it <- read.tree(file.path(MSA, paste0(og, "_tree.txt"))); rt <- read.tree(file.path(RES, paste0(og, "_tree.txt")))
  tr <- root(it, outgroup = clade(rt, rt$edge[rt$edge[, 1] == Ntip(rt) + 1, 2][1]), resolve.root = TRUE, edgelabel = TRUE)
  sp <- sp_of(tr$tip.label); gene <- mapply(function(t, s) norm(s, substring(t, nchar(s) + 2)), tr$tip.label, sp)
  mem <- tr$tip.label[unlist(lapply(names(STUDY), function(x) which(sp == STUDY[[x]] & gene == strip(cl[[x]]))))]
  col <- ifelse(tr$tip.label %in% mem, "#D55E00", ifelse(sp %in% ANGIO, "#0072B2", ifelse(sp %in% GYMNO, "#009E73", "grey40")))
  pdf(file.path(OUT, paste0(og, "_rooted_tree.pdf")), width = 8.5, height = max(6, 0.17 * Ntip(tr)))
  plot(tr, tip.color = col, cex = 0.55, no.margin = TRUE, label.offset = 0.01)
  nodelabels(tr$node.label, frame = "none", cex = 0.45, adj = c(1.1, -0.3))
  legend("topleft", c("angiosperm", "gymnosperm", "moss", "clique member"), text.col = c("#0072B2", "#009E73", "grey40", "#D55E00"),
         bty = "n", cex = 0.7, title = sprintf("%s: IQ-TREE ML tree, UFBoot at nodes, rooted as OrthoFinder roots it", og))
  dev.off()
}
