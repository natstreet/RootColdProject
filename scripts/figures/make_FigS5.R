#!/usr/bin/env Rscript
# make_FigS5.R — Fig. S5: the core's manuscript clusters 1, 2 and 4 (representative spruce and pine members) in the companion conifer
# experiment (van Zalen et al.): root and needle, 5 C, -5 C and drought (ledger N159, N160, R045). Written by SCI; committed by CODE
# with the cluster labels read from Table S3 (as make_Fig5.R), R's pdf device and plotmath degree signs (no reliance on the locale).
# Input: results/a2_per_level_medians.tsv (scripts/companion/a2_tissue_stress_specificity.R), tables/TableS3_cliqueClusterGO_classic.tsv,
#   scripts/figures/FigS5_cluster_labels.tsv (clusters 1 and 2 as Fig. 5; cluster 4, D088). Output: Figure_S5.pdf, Figure_S5.png (300 dpi)
# Usage (working tree): Rscript scripts/figures/make_FigS5.R [input] [out dir]
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE); inp <- if (length(args) > 0) args[1] else "results/a2_per_level_medians.tsv"; outdir <- if (length(args) > 1) args[2] else "."
lab <- fread("scripts/figures/FigS5_cluster_labels.tsv"); t3 <- fread("tables/TableS3_cliqueClusterGO_classic.tsv")
lab[, term := vapply(seq_len(.N), function(i) { h <- t3[`Cluster (manuscript)` == cluster[i] & GO.ID == go_id[i], unique(Term)]
  if (length(h) != 1) stop("GO term ", go_id[i], " is not enriched for cluster ", cluster[i], " in Table S3"); h }, "")]
source("scripts/figures/cluster_labels.R")                                       # approved names (D093)
clab <- sprintf("%s (Cluster %d)", mapply(label_for, lab$cluster, lab$go_id, lab$term), lab$cluster)
d <- fread(inp)
lv <- list("5C" = paste0("5C_", c("6h","24h","3d","10d")), "-5C" = paste0("neg5C_", c("6h","24h","3d","10d")))
lab_cold <- c("6 h","24 h","3 d","10 d"); lab_dr <- c("60%","40%","30%","30%, 7 d","Collapse","Collapse, 2 d")
d[, step := mapply(function(t, l) if (t == "drought") match(l, c("60%","40%","30%","30%7d","Collapse","C2d","FC60","FC40","FC30","FC30d7","Collapsed","Collapsed2d")) %% 6 else match(l, lv[[t]]), trt, level)]
d[trt == "drought" & step == 0, step := 6]
stopifnot(!anyNA(d$step), nrow(d) == 2 * 2 * (4 + 4 + 6))
mc <- paste0("c", lab$cluster, "_median"); stopifnot(all(mc %in% names(d)))
ctl <- unique(d[, .(sp, tissue, trt)])[, level := "control"][, step := 0]
for (k in mc) set(ctl, j = k, value = 0)
d <- rbind(d, ctl, use.names = TRUE, fill = TRUE)
m <- melt(d, id.vars = c("sp","tissue","trt","level","step"), measure.vars = mc, variable.name = "cl", value.name = "lfc")
m[, cl := factor(cl, mc, clab)]
m[, trt := factor(trt, c("5C","-5C","drought"), c("5*degree*C", "-5*degree*C", "Drought"))]   # parsed by label_parsed
m[, row := factor(paste(ifelse(sp == "spruce", "Norway spruce", "Scots pine"), tissue), c("Norway spruce root","Norway spruce needle","Scots pine root","Scots pine needle"))]
labs_all <- c("C", "6 h", "24 h", "3 d", "10 d", "60%", "40%", "30%", "30%, 7 d", "Collapse", "Collapse, 2 d")
m[, xl := factor(ifelse(step == 0, "C", ifelse(trt == "Drought", lab_dr[pmax(step, 1)], lab_cold[pmax(step, 1)])), levels = labs_all)]
p <- ggplot(m, aes(xl, lfc, colour = cl, shape = cl, group = cl)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
  facet_grid(row ~ trt, scales = "free_x", space = "free_x", labeller = labeller(trt = label_parsed)) +
  scale_colour_manual(values = c("#000000", "#999999", "#D55E00")[seq_along(mc)]) +                # Okabe-Ito vermillion: cluster 4
  labs(x = "Time (cold) or soil water content as % of field capacity (drought); C, control",
       y = expression("Median " * log[2] * " fold change"), colour = NULL, shape = NULL) +
  theme_classic(base_size = 8) + theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(size = 7),
                                       axis.text.x = element_text(angle = 45, hjust = 1, size = 6), panel.border = element_rect(fill = NA, colour = "grey80"))
ggsave(file.path(outdir, "Figure_S5.pdf"), p, width = 180, height = 160, units = "mm", device = grDevices::pdf, useDingbats = FALSE)
ggsave(file.path(outdir, "Figure_S5.png"), p, width = 180, height = 160, units = "mm", dpi = 300, bg = "white")
