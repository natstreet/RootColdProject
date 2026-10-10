#!/usr/bin/env Rscript
# make_FigS1.R — Figure S1: worked example of one conserved-core clique (orthogroup given as the first argument; the
# authors' choice of example) and its conserved co-expression neighbourhood.
#   (a) Cold response of the clique member in each genotype, as log2 fold change against the matched control; an open
#       ring marks a time point at which the member is differentially expressed (adjusted P <= 0.05, |log2FC| >= 1).
#   (b) Per genotype, the expression over the time course (median per condition, z-scored per gene) of the member's
#       conserved neighbours: genes in the overlap of its ComPlEx neighbourhood with that of the other member in at
#       least MIN_PAIRS of its five pairwise comparisons (scripts/analysis/core_module_directions.R). Neighbours are
#       coloured by the sign of their correlation with the member over all samples (ComPlEx uses |r|, so both occur);
#       the member is the bold black line and the dashed lines are the mean profiles of the two groups.
#   (c) GO enrichment (biological process) of the Col-0 member's conserved neighbours, with the Arabidopsis
#       annotation: topGO classic Fisher test, Benjamini-Hochberg adjustment, padj < 0.05, as in Tables S1-S3. The
#       background is the Col-0 genes that could be conserved neighbours: genes passing ComPlEx's expression filters
#       with an orthologue in at least MIN_PAIRS of the five other genotypes (scripts/ComPlEx/ComPlEx.R).
#
# Inputs (repository root):
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS   the conserved-core cliques
#   DE/de_results.tsv                         fold changes and DEG calls (scripts/DE/differential_expression.R)
#   tables/core_conserved_neighbours.tsv.gz   conserved neighbours and their correlation with the member
#   results/orthology_class.tsv               the orthogroup's orthology class (scripts/analysis/orthology_class.R, D043)
#                                             (scripts/analysis/core_module_directions.R; its MIN_PAIRS is used here)
#   data/expression/med_expression/<sp>_med.txt.gz         median expression per condition
#   data/expression/<sp>_root_cold_stress_expression.txt   ComPlEx input matrices (GO background)
#   data/ComPlEx/orthologs-<sp1>-<sp2>.RData                ComPlEx orthologue pairs (GO background)
#   data/annotation/Athal_go_ids.tsv, data/annotation/gene_aliases_20140331.txt
# Output: Figure_S1.pdf, Figure_S1.png (300 dpi), Figure_S1_values.tsv (plotted values), Figure_S1_GO.tsv (all
#   significant terms)
# Usage:  Rscript scripts/figures/make_FigS1.R <orthogroup>
#
# ORTHOGROUP IDENTIFIERS ARE RUN-SPECIFIC. They refer to the orthology in data/orthologs/ and are not comparable
# with identifiers from any other OrthoFinder run.

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table) })
gzread <- function(f) fread(cmd = paste("gzip -dc", shQuote(f)))   # gzipped tables without R.utils
source("scripts/GO/go_table_utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript scripts/figures/make_FigS1.R <orthogroup>")
OG <- args[1]
SP <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
LAB <- c(col0 = "Col-0", ost0 = "Ost-0", aspen = "aspen", birch = "birch", spruce = "spruce", pine = "pine")
source("scripts/figures/genotype_colours.R")
COL <- GENOTYPE_COLOURS                                                  # the shared genotype palette
TPS <- c("6h", "24h", "3d", "10d"); TPLAB <- c("6 h", "24 h", "3 d", "10 d")
CONDLAB <- c("control", TPLAB)
POS <- "#8C8C8C"; NEG <- "#E69F00"; GREY <- "#666666"     # grey / Okabe-Ito orange: distinct from the genotype palette
N_TERMS <- 10                                                             # GO terms drawn in panel c
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
strip <- function(x) sub("[.][0-9]+$", "", as.character(x))

core <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
row <- core[core$OrthoGroup == OG, ]
if (nrow(row) != 1) stop(OG, " is not a conserved-core orthogroup")
mem <- unlist(row[SP])
al <- read.delim("data/annotation/gene_aliases_20140331.txt", colClasses = "character")
sym <- al$symbol[toupper(al$locus_name) == strip(mem[["col0"]])]
has_symbol <- length(sym) > 0 && nzchar(sym[1])
label <- if (has_symbol) sym[1] else strip(mem[["col0"]])
ocl <- read.delim("results/orthology_class.tsv")
ocl <- ocl$orthology_class[ocl$OrthoGroup == OG]
if (length(ocl) != 1) stop(OG, " has no orthology class")

## ---- (a) fold changes of the members ----
de <- fread("DE/de_results.tsv")
pa <- rbindlist(lapply(SP, function(s) {
  r <- de[species == s & gene == mem[[s]]][match(TPS, timepoint)]
  if (anyNA(r$gene)) stop("fold changes missing for ", s, " ", mem[[s]])
  data.table(species = s, gene = mem[[s]], timepoint = TPS, x = seq_along(TPS), lfc = r$log2FoldChange, deg = r$deg)
}))
pa[, species := factor(species, SP)]

## ---- (b) conserved neighbours ----
cn <- gzread("tables/core_conserved_neighbours.tsv.gz")[OrthoGroup == OG]
MIN_PAIRS <- unique(cn$min_pairs); stopifnot(length(MIN_PAIRS) == 1)
if (!setequal(cn$species, SP)) stop("conserved neighbours missing for ", paste(setdiff(SP, cn$species), collapse = ", "))
z <- function(m) { m <- t(scale(t(m))); m[is.na(m)] <- 0; m }
pb <- rbindlist(lapply(SP, function(s) {
  x <- gzread(sprintf("data/expression/med_expression/%s_med.txt.gz", s))
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  if (ncol(m) != 5) stop(s, ": expected five conditions in the median expression table")
  c_s <- cn[species == s]
  g <- c(mem[[s]], c_s$neighbour); if (!all(g %in% rownames(m))) stop(s, ": genes missing from median expression")
  zm <- z(m[g, , drop = FALSE])
  data.table(species = s, gene = rep(g, each = 5), x = rep(1:5, length(g)), value = as.vector(t(zm)),
             role = rep(c("member", ifelse(c_s$r > 0, "positive", "negative")), each = 5))
}))
pb[, species := factor(species, SP)]
cnt <- cn[, .(pos = sum(r > 0), neg = sum(r < 0)), by = species]
strip_lab <- setNames(sprintf("%s: %d positive, %d negative", LAB[cnt$species], cnt$pos, cnt$neg), cnt$species)
mean_lines <- pb[role != "member", .(value = mean(value)), by = .(species, role, x)]

## ---- (c) GO enrichment of the Col-0 member's conserved neighbours ----
expressed <- function(sp) {                                   # ComPlEx.R's expression filters
  x <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp))
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  m <- m[matrixStats::rowSds(m) > 0 & matrixStats::rowCounts(m != 0) >= 3, , drop = FALSE]
  rownames(m)[rowMeans(m) > 1.0]
}
e_col0 <- expressed("col0")
nb_ortho <- table(unlist(lapply(setdiff(SP, "col0"), function(t) {            # orthologue in t, both expressed
  o <- ld(sprintf("data/ComPlEx/orthologs-col0-%s.RData", t)); e_t <- expressed(t)
  unique(o$Species1[o$Species1 %in% e_col0 & o$Species2 %in% e_t]) })))
bg <- names(nb_ortho)[nb_ortho >= MIN_PAIRS]
study <- cn[species == "col0", neighbour]
if (!all(study %in% bg)) stop("conserved neighbours outside the GO background")
g2go <- read_go_mapping("data/annotation/Athal_go_ids.tsv")
g2 <- g2go[names(g2go) %in% strip(bg)]
gl <- factor(as.integer(names(g2) %in% strip(study)), levels = 0:1); names(gl) <- names(g2)
GOdata <- suppressMessages(new("topGOdata", ontology = "BP", allGenes = gl, annot = annFUN.gene2GO, gene2GO = g2))
res <- suppressMessages(runTest(GOdata, algorithm = "classic", statistic = "fisher"))
sig <- go_result(GOdata, res)
fwrite(sig, "Figure_S1_GO.tsv", sep = "\t")
pc <- head(sig, N_TERMS)
pc$Term <- factor(pc$Term, rev(pc$Term))

## ---- drawing ----
th <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 6, colour = "black"), plot.title = element_text(size = 8, hjust = 0.5),
        axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        plot.tag = element_text(face = "bold", size = 10))

A <- ggplot(pa, aes(x, lfc, colour = species)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.4) +
  geom_line(linewidth = 0.8, na.rm = TRUE) + geom_point(size = 1.6, na.rm = TRUE) +
  geom_point(data = pa[deg & !is.na(lfc)], shape = 21, size = 4, stroke = 0.7, fill = NA, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.2, size = 2.1, colour = GREY,
           label = "'ring: DE (adjusted'~italic(P)<=0.05*','~abs(log[2]*FC)>=1*')'", parse = TRUE) +
  scale_colour_manual(values = COL, labels = LAB, name = NULL) +
  scale_x_continuous(breaks = seq_along(TPS), labels = TPLAB, expand = expansion(mult = 0.06)) +
  scale_y_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), 1), expand = expansion(mult = 0.1)) +
  labs(x = NULL, y = expression(log[2]~fold~change~vs~control), tag = "a",
       title = if (has_symbol)                      # gene symbol italic, locus identifier roman (NRS, ledger R049)
         bquote(atop(italic(.(label)) ~ "clique members (" * .(OG) * ")", .(paste(ocl, "(gene-tree call)"))))
       else sprintf("%s clique members (%s)\n%s (gene-tree call)", label, OG, ocl)) + th +
  theme(legend.position = "bottom", legend.text = element_text(size = 7), legend.key.width = unit(10, "pt"),
        legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8)))

B <- ggplot(pb[role != "member"], aes(x, value, group = gene, colour = role)) +
  geom_line(linewidth = 0.2, alpha = 0.25) +
  geom_line(data = mean_lines, aes(x, value, group = role, colour = role), linewidth = 0.8, linetype = "22",
            inherit.aes = FALSE) +
  geom_line(data = pb[role == "member"], colour = "black", linewidth = 1.1) +
  facet_wrap(~ species, nrow = 2, labeller = labeller(species = strip_lab)) +
  scale_colour_manual(values = c(positive = POS, negative = NEG), name = NULL,
                      labels = c(positive = "positively correlated with the member",
                                 negative = "negatively correlated with the member")) +
  scale_x_continuous(breaks = 1:5, labels = CONDLAB, expand = expansion(mult = 0.04)) +
  labs(x = NULL, y = "expression (z-score of condition medians)", tag = "b",
       title = sprintf("Conserved neighbours of each member (in the overlap in at least %d of 5 comparisons); black: member; dashed: group means",
                       MIN_PAIRS)) +
  th + theme(strip.background = element_blank(), strip.text = element_text(size = 7), panel.spacing.x = unit(14, "pt"),
             legend.position = "bottom", legend.text = element_text(size = 7))

C <- ggplot(pc, aes(-log10(padj), Term)) +
  geom_col(fill = COL[["col0"]], width = 0.65) +
  geom_text(aes(label = sprintf("%d/%d", Significant, Annotated)), hjust = -0.1, size = 2.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression(-log[10]~adjusted~italic(P)), y = NULL, tag = "c",
       title = sprintf("GO (BP) of the Col-0 member's %d conserved neighbours\nlabels: neighbours / background genes with the term",
                       length(study))) + th + theme(plot.title.position = "plot")

fig <- (A | C) / B + plot_layout(heights = c(0.85, 1.25))
ggsave("Figure_S1.pdf", fig, width = 7.2, height = 7.4, bg = "white")
ggsave("Figure_S1.png", fig, width = 7.2, height = 7.4, dpi = 300, bg = "white")

vals <- rbind(
  data.table(panel = "a", species = as.character(pa$species), quantity = paste0("log2FC_", pa$timepoint),
             value = pa$lfc, note = ifelse(pa$deg, "DE", "")),
  data.table(panel = "b", species = cnt$species, quantity = "conserved_neighbours_positive", value = cnt$pos, note = ""),
  data.table(panel = "b", species = cnt$species, quantity = "conserved_neighbours_negative", value = cnt$neg, note = ""),
  data.table(panel = "c", species = "col0", quantity = paste("padj", pc$GO.ID), value = pc$padj, note = as.character(pc$Term)))
fwrite(vals, "Figure_S1_values.tsv", sep = "\t")
cat(sprintf("wrote Figure_S1.pdf/.png for %s (%s); conserved neighbours in >= %d of 5 comparisons\n", OG, label, MIN_PAIRS))
print(cnt)
cat(sprintf("GO: %d of %d background genes are study genes with GO; %d significant terms; top: %s\n",
            sum(gl == 1), length(gl), nrow(sig), paste(head(sig$Term, 3), collapse = "; ")))
