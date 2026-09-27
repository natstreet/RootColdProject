#!/usr/bin/env Rscript
# make_FigS6.R — Figure S6: robustness of the conserved co-expression core.
#   (a) Conservation permutation test: the observed number of orthogroups with a co-expressolog in all 15 comparisons
#       against the label-permutation null.
#   (b) Cluster coherence in independent data (STRING): observed within-cluster co-expressed gene pairs among the core
#       Arabidopsis genes against the label-permuted null.
#
# This step only plots. The null distributions and observed values are produced by:
#   panel (a): scripts/validation/conservation_label_permutation.R -> conservation_null_null.csv,
#              conservation_null_obs.txt   (from the figshare deposit)
#   panel (b): scripts/validation/make_core_clusters.R, then string_external_validation.R -> string_null.csv,
#              string_obs.txt   (needs the STRING v12 A. thaliana download)
# Run those first, then run this from the same directory.
# Output: Figure_S6.pdf, Figure_S6.png (300 dpi)
# Usage:  Rscript scripts/figures/make_FigS6.R

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

nullA <- read.csv("conservation_null_null.csv")$null_core_size; obsA <- as.integer(readLines("conservation_null_obs.txt"))
nullB <- read.csv("string_null.csv")$null;                     obsB <- as.integer(readLines("string_obs.txt"))
pA <- (sum(nullA >= obsA) + 1) / (length(nullA) + 1)
pB <- (sum(nullB >= obsB) + 1) / (length(nullB) + 1)
foldB <- obsB / max(mean(nullB), 1e-9)
RED <- "#d62728"

th <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text = element_text(colour = "black"),
        plot.tag = element_text(face = "bold", size = 14))

A <- ggplot(data.frame(x = nullA), aes(x)) +
  geom_histogram(binwidth = 1, boundary = -0.5, fill = "#add8e6") +
  geom_vline(xintercept = obsA, colour = RED, linewidth = 1.2) +
  annotate("segment", x = obsA * 0.6, y = length(nullA) * 0.6, xend = obsA - 1, yend = 5, colour = RED,
           arrow = arrow(length = unit(6, "pt"))) +
  annotate("text", x = obsA * 0.55, y = length(nullA) * 0.62, label = sprintf("observed = %d", obsA), colour = RED,
           size = 4.6, hjust = 0.5, vjust = -0.3) +
  coord_cartesian(xlim = c(-3, obsA + 6)) + scale_x_continuous(breaks = seq(0, obsA + 6, 50)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Orthogroups conserved across all 15\ncomparisons (permuted null)", y = "Permutations", tag = "a",
       title = sprintf("Conservation permutation test (null max = %d; P = %.4g)", max(nullA), pA)) + th

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

fig <- A | B
ggsave("Figure_S6.pdf", fig, width = 12.5, height = 4.6)
ggsave("Figure_S6.png", fig, width = 12.5, height = 4.6, dpi = 300)
cat(sprintf("wrote Figure_S6.pdf/.png ; panel a P=%.4g, panel b P=%.4g, fold=%.1f\n", pA, pB, foldB))
