#!/usr/bin/env Rscript
# make_FigS4.R — Figure S4: robustness of the conserved co-expression core.
#   (a) Unrelated-partner null (D097, R071, R073): the observed number of orthogroups (i) with co-expressologs in all 15
#       comparisons and (ii) with a complete six-genotype clique after the size rule, against size-matched pseudo-families of genes
#       from other orthogroups. Each sub-panel has its own linear x axis from 0 to beyond the observed value, so the null and the
#       observed value are both visible.
#   (b) Cluster coherence in independent data (STRING): observed within-cluster co-expressed gene pairs among the core
#       Arabidopsis genes against the label-permuted null.
#
# This step only plots. The null distributions and observed values are produced by:
#   panel (a): scripts/validation/pseudo_families.R, pseudo_family_complex.R (micro; 100 replicates x 15 comparisons) and
#              pseudo_family_summary.R -> data/pseudo_family_null/per_replicate.tsv, null_summary.tsv
#   panel (b): scripts/validation/make_core_clusters.R, then string_external_validation.R -> string_null.csv,
#              string_obs.txt   (needs the STRING v12 A. thaliana download)
# Run those first, then run this from the same directory.
# Output: Figure_S4.pdf, Figure_S4.png (300 dpi)
# Usage:  Rscript scripts/figures/make_FigS4.R

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

pr <- read.delim("data/pseudo_family_null/per_replicate.tsv"); ns <- read.delim("data/pseudo_family_null/null_summary.tsv")
tot <- aggregate(n ~ replicate + statistic, pr, sum)
nullB <- read.csv("string_null.csv")$null;                     obsB <- as.integer(readLines("string_obs.txt"))
pB <- (sum(nullB >= obsB) + 1) / (length(nullB) + 1)
foldB <- obsB / max(mean(nullB), 1e-9)
RED <- "#d62728"

th <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text = element_text(colour = "black"),
        plot.tag = element_text(face = "bold", size = 14))

panelA <- function(stat, xlab, tag) {
  v <- tot$n[tot$statistic == stat]; o <- ns$observed[ns$statistic == stat]; nrep <- ns$replicates[ns$statistic == stat]
  stopifnot(length(v) == nrep, length(o) == 1)
  ggplot(data.frame(x = v), aes(x)) +
    geom_histogram(binwidth = 1, boundary = -0.5, fill = "#add8e6") +
    geom_vline(xintercept = o, colour = RED, linewidth = 1.2) +
    annotate("text", x = o, y = Inf, label = sprintf("observed = %d", o), colour = RED, size = 4, hjust = 1.05, vjust = 1.5) +
    coord_cartesian(xlim = c(0, o * 1.05)) + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = xlab, y = "Replicates", tag = tag,
         title = sprintf("Null mean %.1f (range %d-%d; %d replicates)", mean(v), min(v), max(v), nrep)) + th
}
A1 <- panelA("all15", "Orthogroups with co-expressologs in all 15 comparisons", "a")
A2 <- panelA("clique_rule", "Orthogroups with a complete clique (after the size rule)", NULL)
hB <- hist(nullB, breaks = seq(min(nullB), max(nullB), length.out = 31), plot = FALSE)   # 30 equal-width bins
B <- ggplot(data.frame(x = nullB), aes(x)) +
  geom_histogram(breaks = hB$breaks, fill = "#90ee90", colour = "white") +
  geom_vline(xintercept = obsB, colour = RED, linewidth = 1.2) +
  annotate("segment", x = obsB * 0.64, y = max(hB$counts) * 0.6, xend = obsB - 2, yend = max(hB$counts) * 0.02,
           colour = RED, arrow = arrow(length = unit(6, "pt"))) +
  annotate("text", x = obsB * 0.60, y = max(hB$counts) * 0.62, label = sprintf("observed = %d", obsB), colour = RED,
           size = 4.6, hjust = 0.5, vjust = -0.3) +
  coord_cartesian(xlim = c(min(nullB) - 20, obsB + 30)) + scale_x_continuous(breaks = seq(0, obsB + 30, 50)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Within-cluster co-expressed gene pairs\nin STRING (clusters randomised)", y = "Permutations", tag = "b",
       title = sprintf("Cluster coherence in independent data (%.1f-fold; P = %.4g)", foldB, pB)) + th

fig <- (A1 / A2) | B
ggsave("Figure_S4.pdf", fig, width = 12.5, height = 6.4, bg = "white")
ggsave("Figure_S4.png", fig, width = 12.5, height = 6.4, dpi = 300, bg = "white")
cat(sprintf("wrote Figure_S4.pdf/.png ; panel a: %s ; panel b P=%.4g, fold=%.1f\n", paste(sprintf("%s observed %d, null mean %.2f", ns$statistic, ns$observed, ns$null_mean), collapse = "; "), pB, foldB))
