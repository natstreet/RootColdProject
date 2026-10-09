#!/usr/bin/env Rscript
# make_Fig5.R — Figure 5: the two growth components of the conserved core (manuscript clusters 1 and 2) respond separably to cold
# (ledger R043; producer written by SCI, committed by CODE with two changes: the panel labels are read from Table S3, and the
# axis limits follow the data instead of fixed values that could silently drop points).
#   a, b   median log2 fold change over the cluster's clique members, per genotype and time point (0 = control)
#   c      sustained response: per genotype, the mean over 24 h, 3 d and 10 d of each cluster's median log2 fold change
# Labels: for each cluster, the GO term named in scripts/figures/Fig5_cluster_labels.tsv, read with its name from Table S3
#   (tables/TableS3_cliqueClusterGO_classic.tsv), where it must be enriched for that cluster; the script stops otherwise.
# Inputs (working tree): results/cluster_profiles/cluster_profiles.tsv (cluster_profiles.R), tables/TableS3_cliqueClusterGO_classic.tsv
#   (make_TableS3.R), scripts/figures/Fig5_cluster_labels.tsv, scripts/figures/genotype_colours.json
# Outputs: figures/Figure5.pdf, figures/Figure5.png (300 dpi, opaque), results/Figure5_values.tsv
# Usage (working tree): Rscript scripts/figures/make_Fig5.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(jsonlite); library(patchwork) })
cp  <- fread("results/cluster_profiles/cluster_profiles.tsv")
col <- unlist(fromJSON("scripts/figures/genotype_colours.json"))
lab <- fread("scripts/figures/Fig5_cluster_labels.tsv"); t3 <- fread("tables/TableS3_cliqueClusterGO_classic.tsv")
lab[, term := vapply(seq_len(.N), function(i) { h <- t3[`Cluster (manuscript)` == cluster[i] & GO.ID == go_id[i], unique(Term)]
  if (length(h) != 1) stop("GO term ", go_id[i], " is not enriched for cluster ", cluster[i], " in Table S3"); h }, "")]
lab[, display := paste0(toupper(substr(term, 1, 1)), substring(term, 2))]
geno <- c(col0 = "Col-0", ost0 = "Ost-0", aspen = "Aspen", birch = "Birch", spruce = "Norway spruce", pine = "Scots pine")
tps  <- c("0", "6h", "24h", "3d", "10d"); tlab <- c("0", "6 h", "24 h", "3 d", "10 d")
d <- cp[ms_cluster %in% lab$cluster, .(cluster = ms_cluster, species, timepoint, n, median_log2fc)]
stopifnot(nrow(d) == nrow(lab) * 6 * 4, setequal(unique(d$species), names(geno)))
d <- rbind(d, unique(d[, .(cluster, species, n)])[, `:=`(timepoint = "0", median_log2fc = 0)], use.names = TRUE)
d[, t := factor(timepoint, levels = tps, labels = tlab)][, g := factor(geno[species], levels = geno)]
pal <- setNames(col[names(geno)], geno)
pad <- function(x, f = 0.08) { r <- range(x); r + c(-1, 1) * f * diff(r) }
ylim_ab <- pad(c(d$median_log2fc, 0))
th <- theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6), legend.text = element_text(size = 7),
  plot.title = element_text(size = 8, hjust = 0), plot.tag = element_text(face = "bold", size = 10), legend.title = element_blank())
line_panel <- function(k) {
  ggplot(d[cluster == k], aes(t, median_log2fc, colour = g, group = g)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_line(linewidth = 0.6) + geom_point(size = 1.2) +
    scale_colour_manual(values = pal, guide = guide_legend(nrow = 1)) + scale_y_continuous(limits = ylim_ab, breaks = seq(-3, 3, 1)) +
    labs(x = expression("Time at 5 " * degree * "C"), y = expression("Median " * log[2] * " fold change"),
         title = sprintf("%s (Cluster %d, %d orthogroups)", lab[cluster == k, display], k, unique(d[cluster == k]$n))) + th
}
k1 <- lab$cluster[1]; k2 <- lab$cluster[2]
pa <- line_panel(k1); pb <- line_panel(k2) + theme(axis.title.y = element_blank())
w <- dcast(d[timepoint %in% c("24h", "3d", "10d")], g ~ cluster, value.var = "median_log2fc", fun.aggregate = mean)
setnames(w, as.character(c(k1, k2)), c("xa", "xb"))
xl <- pad(c(w$xa, 0), 0.30); yl <- pad(c(w$xb, 0), 0.20)   # room for the point labels
n1 <- tolower(lab[cluster == k1, term]); n2 <- tolower(lab[cluster == k2, term])
pc <- ggplot(w, aes(xa, xb, colour = g)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 2.2) +
  geom_text(data = w[g != "Ost-0"], aes(label = g), size = 2.4, vjust = -1, show.legend = FALSE) +
  geom_text(data = w[g == "Ost-0"], aes(label = g), size = 2.4, hjust = 1.25, show.legend = FALSE) +
  annotate("text", x = xl[1], y = yl[2], label = sprintf("%s down,\n%s up", n1, n2), hjust = 0, vjust = 1, size = 2.3, colour = "grey35") +
  annotate("text", x = xl[2], y = yl[2], label = "both up", hjust = 1, vjust = 1, size = 2.3, colour = "grey35") +
  annotate("text", x = xl[1], y = yl[1], label = "both down", hjust = 0, vjust = 0, size = 2.3, colour = "grey35") +
  scale_colour_manual(values = pal, guide = "none") + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl) +
  labs(x = bquote(.(lab[cluster == k1, display]) ~ "(mean" ~ log[2] ~ "FC)"), y = bquote(.(lab[cluster == k2, display]) ~ "(mean" ~ log[2] ~ "FC)"),
       title = "Sustained response (24 h to 10 d)") + th
fig <- (pa + pb + pc) + plot_layout(widths = c(1, 1, 1.05), guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
dir.create("figures", showWarnings = FALSE); dir.create("results", showWarnings = FALSE)
ggsave("figures/Figure5.pdf", fig, device = grDevices::pdf, width = 180, height = 72, units = "mm", useDingbats = FALSE)   # cairo_pdf needs X11 libraries not present on every system
ggsave("figures/Figure5.png", fig, width = 180, height = 72, units = "mm", dpi = 300, bg = "white")
fwrite(d[order(cluster, g, t), .(cluster, genotype = species, timepoint, n_orthogroups = n, median_log2fc)], "results/Figure5_values.tsv", sep = "\t")
fwrite(w[, .(genotype = g, mean_24h_to_10d_cluster_a = xa, mean_24h_to_10d_cluster_b = xb)], "results/Figure5c_values.tsv", sep = "\t")
cat("Figure 5:", nrow(d), "values; panel labels:", paste(lab$display, collapse = ", "), "\n")
