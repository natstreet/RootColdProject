#!/usr/bin/env Rscript
# motif_conservation_test.R — do co-expressolog pairs share more promoter motifs than non-co-expressolog pairs of the
# same orthogroup? (DECISIONS.md D062 item 4, D063; design R040 as amended in A035.)
#
# Criteria (fixed before the run, A035):
#   motif set   per gene and window: the JASPAR 2024 plant CORE clusters (jaspar_clusters.py) with at least one FIMO hit
#               (p < 1e-4, motif_job_micro.sh); scored per cluster so a family's near-identical matrices count once
#   measure     per pair: Jaccard similarity of the two genes' cluster sets (0 if both are empty)
#   per OG      in each comparison: mean Jaccard of co-expressolog pairs minus mean of non-co-expressolog pairs, when the
#               orthogroup has both pair types there (pairs as in motif_pairs.R)
#   primary     the 6 angiosperm-gymnosperm comparisons (col0, aspen, birch x spruce, pine): per orthogroup, the mean of its
#               per-comparison differences; Wilcoxon signed-rank across orthogroups (two-sided, exact = FALSE)
#   outcome     "supported" if the Hodges-Lehmann estimate is positive and P < 0.05; "not supported" otherwise; the 95% CI of the
#               estimate (Wilcoxon) is reported so a null can be read as the size of difference it excludes
#   secondary   the same per comparison for all 10 comparisons (BH across the 10), shown with the 4 within-lineage ones
#   start codon (D064) genes whose model has no annotated start codon (check_start_codons.py, start_codon FALSE) are
#               dropped with their pairs in mode "primary"; mode "sensitivity" keeps them
# Orthogroups: those listed in <core list> (column OrthoGroup), i.e. the final core after R039.
# Inputs: <pairs> (motif_pairs.tsv.gz), <hits dir> (hits_<genotype>_<win>.tsv.gz, <genotype>_<win>_stats.tsv,
#   start_codons.tsv), <clusters tsv>, <core list>, <win>, <mode> (primary | sensitivity)
# Outputs (<out dir>): motif_pair_jaccard_<win>_<mode>.tsv.gz, motif_og_differences_<win>_<mode>.tsv, motif_test_<win>_<mode>.tsv
# Usage: Rscript motif_conservation_test.R <pairs> <hits dir> <clusters tsv> <core list> <win> <mode> <out dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
PAIRS <- a[1]; HITS <- a[2]; CL <- a[3]; CORE <- a[4]; WIN <- as.integer(a[5]); MODE <- a[6]; OUT <- a[7]
if (!MODE %in% c("primary", "sensitivity")) stop("mode must be primary or sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
G <- c("col0", "aspen", "birch", "spruce", "pine"); ANG <- c("col0", "aspen", "birch"); GYM <- c("spruce", "pine")
cl <- fread(CL)
core <- unique(fread(CORE)$OrthoGroup)
pr <- fread(cmd = paste("gzip -dc", shQuote(PAIRS)))[OrthoGroup %in% core]
sc <- fread(file.path(HITS, "start_codons.tsv"))
nost <- sc[start_codon == FALSE, paste(genotype, gene)]
pr[, no_start := paste(s1, gene1) %in% nost | paste(s2, gene2) %in% nost]
cat("pairs with a gene lacking an annotated start codon:", pr[no_start == TRUE, .N], "of", nrow(pr), "\n")
if (MODE == "primary") pr <- pr[no_start == FALSE]
sets <- rbindlist(lapply(G, function(g) {
  h <- fread(cmd = paste("gzip -dc", shQuote(file.path(HITS, sprintf("hits_%s_%d.tsv.gz", g, WIN)))))
  if (h[!matrix_id %in% cl$matrix_id, .N]) stop("matrices without cluster in ", g)
  unique(merge(h, cl[, .(matrix_id, cluster)], by = "matrix_id")[, .(genotype = g, gene, cluster)]) }))
genes_in <- function(g) unique(c(pr[s1 == g, gene1], pr[s2 == g, gene2]))
for (g in G) { m <- setdiff(genes_in(g), fread(file.path(HITS, sprintf("%s_%d_stats.tsv", g, WIN)))$gene)
  if (length(m)) stop(length(m), " ", g, " genes have no promoter sequence, e.g. ", m[1]) }
S <- split(sets$cluster, paste(sets$genotype, sets$gene))
jac <- function(x, y) { u <- length(union(x, y)); if (u == 0) 0 else length(intersect(x, y)) / u }
pr[, jaccard := mapply(function(a1, g1, a2, g2) jac(S[[paste(a1, g1)]], S[[paste(a2, g2)]]), s1, gene1, s2, gene2)]
fwrite(pr, file.path(OUT, sprintf("motif_pair_jaccard_%d_%s.tsv.gz", WIN, MODE)), sep = "\t")
d <- pr[, .(n_co = sum(coexpressolog), n_non = sum(!coexpressolog), diff = mean(jaccard[coexpressolog]) - mean(jaccard[!coexpressolog])),
        by = .(s1, s2, OrthoGroup)][n_co > 0 & n_non > 0]
d[, angio_gymno := (s1 %in% ANG & s2 %in% GYM) | (s2 %in% ANG & s1 %in% GYM)]
fwrite(d, file.path(OUT, sprintf("motif_og_differences_%d_%s.tsv", WIN, MODE)), sep = "\t")
wt <- function(x) { if (length(x) < 2) return(list(n = length(x), estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
  t <- wilcox.test(x, mu = 0, exact = FALSE, conf.int = TRUE)
  list(n = length(x), estimate = unname(t$estimate), ci_low = t$conf.int[1], ci_high = t$conf.int[2], p = t$p.value) }
pool <- d[angio_gymno == TRUE, .(diff = mean(diff), comparisons = .N), by = OrthoGroup]
pt <- wt(pool$diff)
res <- rbind(data.table(test = "primary: angiosperm-gymnosperm pooled", comparisons = 6L, orthogroups = pt$n, hl_estimate = pt$estimate,
                        ci_low = pt$ci_low, ci_high = pt$ci_high, p = pt$p, p_bh = NA_real_,
                        outcome = ifelse(!is.na(pt$p) & pt$estimate > 0 & pt$p < 0.05, "supported", "not supported")),
             d[, { t <- wt(diff); .(test = paste("secondary:", s1, s2), comparisons = 1L, orthogroups = t$n, hl_estimate = t$estimate,
                                   ci_low = t$ci_low, ci_high = t$ci_high, p = t$p) }, by = .(s1, s2)][
               , p_bh := p.adjust(p, "BH")][, outcome := NA_character_][, !c("s1", "s2")], use.names = TRUE)
res[, `:=`(window = WIN, mode = MODE)]
fwrite(res, file.path(OUT, sprintf("motif_test_%d_%s.tsv", WIN, MODE)), sep = "\t"); print(res, row.names = FALSE)
