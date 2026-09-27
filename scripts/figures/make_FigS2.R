#!/usr/bin/env Rscript
# make_FigS2.R — Figure S2: worked example of resource use, conserved-core orthogroup OG0000118 (NAC transcription
# factors, ATAF/NAM clade).
#   (a) Cold response of the orthogroup's member gene in each genotype across the time course, as log2 fold change
#       against the matched control; an open ring marks a timepoint at which it is differentially expressed (adjusted
#       P <= 0.05 and |log2FC| >= 1). The panel title is built from those DEG calls.
#   (b) The six members and their co-expressolog relationships between species.
#   (c) Number of genes co-expressed with each member in that species' whole-genome co-expression network, with the
#       transcription-factor subset overlaid.
#
# Inputs, all from the figshare deposit (version 4 or later; the open-format tables and the per-species networks were
# introduced in version 4). Run from the directory that contains data/ :
#   data/open/conserved_core.tsv          the 116 conserved-core orthogroups
#   data/open/gene_index.tsv              per-gene fold changes
#   data/open/co_expressologs.tsv         co-expressolog pairs, used for panel (b)
#   data/open/transcription_factors.tsv   TF assignments, used for panel (c)
#   data/networks/network_<species>_hrr100.tsv   per-species co-expression networks
#   DE/<species>/DEGs_<species>_C_vs_<timepoint>.txt   per-timepoint DEG lists written by
#                                         scripts/DE/differential_expression.R (the rings in panel a)
# Output: Figure_S2.pdf, Figure_S2.png (300 dpi), Figure_S2_values.tsv (the plotted values)
#
# ORTHOGROUP IDENTIFIERS ARE RUN-SPECIFIC. OG0000118 refers to the orthology in data/orthologs/; it is not comparable
# with identifiers from any other OrthoFinder run.

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

OG <- "OG0000118"
SP <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
LAB <- c(col0 = "Col-0", ost0 = "Ost-0", aspen = "aspen", birch = "birch", spruce = "spruce", pine = "pine")
COL <- c(col0 = "#0072B2", ost0 = "#56B4E9", aspen = "#009E73", birch = "#7FBF3F", spruce = "#D55E00", pine = "#CC79A7")
TPS <- c("6h", "24h", "3d", "10d"); TPLAB <- c("6 h", "24 h", "3 d", "10 d")
GREY <- "#666666"
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, quote = "", check.names = FALSE)

core <- rd("data/open/conserved_core.tsv"); gi <- rd("data/open/gene_index.tsv")
coex <- rd("data/open/co_expressologs.tsv"); tf <- rd("data/open/transcription_factors.tsv")
row <- core[core$OrthoGroup == OG, ]; stopifnot(nrow(row) == 1)
mem <- unlist(row[SP]); mem <- mem[!is.na(mem) & nzchar(mem)]
tfset <- paste(tf$species, tf$gene)

# Differentially expressed at a timepoint as defined in the Methods (Wald test against no change, adjusted P <= 0.05
# and |log2FC| >= 1): membership of that timepoint's DEG list from differential_expression.R.
is_deg <- function(sp, gene, tp) {
  f <- list.files(file.path("DE", sp), pattern = sprintf("^DEGs_%s_C_vs_%s", sp, tp), full.names = TRUE)
  if (length(f) != 1) stop("expected one DEG list for ", sp, " ", tp, " in DE/", sp)
  gene %in% scan(f, what = "", quiet = TRUE)
}

# panel (a) values
pa <- do.call(rbind, lapply(names(mem), function(s) {
  r <- gi[gi$species == s & gi$gene == mem[[s]], ]; stopifnot(nrow(r) == 1)
  data.frame(species = s, gene = mem[[s]], timepoint = TPS, x = seq_along(TPS),
             lfc = unlist(r[paste0("lfc_", TPS)]), deg = vapply(TPS, function(t) is_deg(s, mem[[s]], t), TRUE))
}))
pa$species <- factor(pa$species, SP)

# title from the DEG calls: each differentially expressed genotype, by the sign of its DEG fold changes
join <- function(x) if (length(x) == 1) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
dirs <- list()
for (s in SP) {
  v <- pa$lfc[pa$species == s & pa$deg & !is.na(pa$lfc)]
  if (length(v)) { k <- if (max(v) < 0) "repressed" else if (min(v) > 0) "induced" else "changed"; dirs[[k]] <- c(dirs[[k]], LAB[[s]]) }
}
ttl <- paste(vapply(names(dirs), function(k) paste(k, "in", join(dirs[[k]])), ""), collapse = "; ")
ttl <- paste0(toupper(substr(ttl, 1, 1)), substring(ttl, 2))

# panel (b): species pairs linked by a co-expressolog between the six members
sel <- coex[coex$OrthologGroup == OG & coex$GeneSpecies1 %in% mem & coex$GeneSpecies2 %in% mem, ]
pairs <- unique(t(apply(sel[, c("Species1", "Species2")], 1, sort)))
ang <- (seq(90, 450, length.out = 7)[-7]) * pi / 180
pos <- data.frame(species = factor(SP, SP), x = cos(ang), y = sin(ang))
edges <- data.frame(x = pos$x[match(pairs[, 1], SP)], y = pos$y[match(pairs[, 1], SP)],
                    xend = pos$x[match(pairs[, 2], SP)], yend = pos$y[match(pairs[, 2], SP)])

# panel (c): network neighbours and the TF subset
pc <- do.call(rbind, lapply(names(mem), function(s) {
  n <- rd(sprintf("data/networks/network_%s_hrr100.tsv", s)); g <- mem[[s]]
  sub <- n[n$source == g | n$target == g, ]
  partner <- ifelse(sub$source == g, sub$target, sub$source)
  data.frame(species = s, co_expressed = nrow(sub), transcription_factors = sum(paste(s, partner) %in% tfset))
}))
pc$species <- factor(pc$species, rev(SP))

th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 6, colour = "black"), plot.title = element_text(size = 8, hjust = 0.5),
        axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        plot.tag = element_text(face = "bold", size = 10))

A <- ggplot(pa, aes(x, lfc, colour = species)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.4) +
  geom_line(linewidth = 0.8, na.rm = TRUE) + geom_point(size = 1.6, na.rm = TRUE) +
  geom_point(data = pa[pa$deg & !is.na(pa$lfc), ], shape = 21, size = 4, stroke = 0.7, fill = NA, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.2, size = 2.1, colour = GREY,
           label = "'open ring: differentially expressed (adjusted'~italic(P)<=0.05*','~abs(log[2]*FC)>=1*')'", parse = TRUE) +
  scale_colour_manual(values = COL, labels = LAB, name = NULL) +
  scale_x_continuous(breaks = seq_along(TPS), labels = TPLAB, expand = expansion(mult = 0.06)) +
  scale_y_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), 1), expand = expansion(mult = 0.1)) +
  labs(x = NULL, y = expression(log[2]~fold~change~vs~control), title = ttl, tag = "a") + th +
  theme(legend.position = "bottom", legend.text = element_text(size = 7), legend.key.width = unit(10, "pt"),
        legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8)))

B <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend), colour = "#BBBBBB", linewidth = 0.45) +
  geom_point(data = pos, aes(x, y, colour = species), size = 6, show.legend = FALSE) +
  geom_text(data = pos, aes(x * 1.5, y * 1.5, label = LAB[as.character(species)]), size = 2.1) +
  scale_colour_manual(values = COL) + coord_equal(xlim = c(-1.85, 1.85), ylim = c(-1.7, 1.7)) +
  labs(title = sprintf("Co-expressologs between the members:\n%d of 15 species pairs linked", nrow(pairs)), tag = "b") +
  theme_void(base_size = 8) + theme(plot.title = element_text(size = 8, hjust = 0.5),
                                    plot.tag = element_text(face = "bold", size = 10))

C <- ggplot(pc, aes(y = species)) +
  geom_col(aes(x = co_expressed, fill = species), width = 0.62, show.legend = FALSE) +
  geom_col(aes(x = transcription_factors), fill = "#333333", width = 0.26) +
  geom_text(aes(x = co_expressed + 2, label = co_expressed), hjust = 0, size = 2.1) +
  scale_fill_manual(values = COL) + scale_y_discrete(labels = LAB) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = "co-expressed genes", y = NULL, tag = "c",
       title = "Genes co-expressed with each member;\ndark = transcription factors") + th

fig <- A / (B | C) + plot_layout(heights = c(1, 0.95))
ggsave("Figure_S2.pdf", fig, width = 7.2, height = 5.4)
ggsave("Figure_S2.png", fig, width = 7.2, height = 5.4, dpi = 300)

vals <- rbind(
  data.frame(panel = "a", species = as.character(pa$species), quantity = paste0("log2FC_", pa$timepoint),
             value = pa$lfc, deg = pa$deg),
  data.frame(panel = "b", species = NA, quantity = "species_pairs_linked", value = nrow(pairs), deg = NA),
  data.frame(panel = "c", species = as.character(pc$species), quantity = "co_expressed_genes", value = pc$co_expressed, deg = NA),
  data.frame(panel = "c", species = as.character(pc$species), quantity = "of_which_TFs", value = pc$transcription_factors, deg = NA))
write.table(vals, "Figure_S2_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("wrote Figure_S2.pdf/.png for %s\n  species pairs linked: %d of 15\n", OG, nrow(pairs)))
cat("  co-expressed genes:", paste(sprintf("%s=%d", pc$species, pc$co_expressed), collapse = " "), "\n")
cat("  of which transcription factors:", paste(sprintf("%s=%d", pc$species, pc$transcription_factors), collapse = " "), "\n")
cat("  title:", ttl, "\n")
