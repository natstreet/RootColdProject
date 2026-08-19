#!/usr/bin/env Rscript
# make_Fig1.R — Figure 1: orthogroup overlap of expressed genes (UpSet, panel a), and the
# composition of expressed genes (b) and differentially expressed genes (c, d) in the
# orthogroups shared across all species/ecotypes.
#
# Reproduces the published figure from the deposited data plus two intermediate inputs that are
# provided in this repository (they are not in the figshare deposit): the per-timepoint MEDIAN
# expression tables and Elena van Zalen's per-species gene->orthogroup maps. The published figure
# was assembled externally in Illustrator; this script regenerates its underlying values and a
# composite rendering.
#
# Definitions (matching the original analysis):
#   * expressed in a condition: median expression > 0 in the control timepoint, or in >=1 of the
#     cold timepoints (6h/24h/3d/10d), per species. An orthogroup is expressed in a condition if
#     >=1 of its genes is. Panel a = UpSet of the 12 species x {control, treatment} sets; the
#     top intersection (all 12) is the shared expressed core.
#   * panel b: expressed genes (in the GO background bg_*) whose orthogroup is in that shared core.
#   * panels c/d: DEGs categorised by how many species have a DEG in the orthogroup — "DEG across
#     all" (all six), "in some" (2-5), "species-specific" (1), "singletons" (no orthogroup). Panel
#     d shows, per species, the DEGs in the "DEG across all" orthogroups.
#
# Inputs (run from a directory with the unpacked figshare deposit `data/`, into which this repo's
# committed data/ files merge):
#   data/expression/med_expression/<sp>_med.txt.gz          per-timepoint median expression (repo)
#   data/ComPlEx/gene2OG/gene2OG_<Sp>.rda                    gene -> orthogroup maps (repo)
#   data/annotation/bg_<sp>_<code>.rda                       expressed-gene backgrounds (repo)
#   data/DEGs/DEGs_<sp>.RData                                per-species DEG lists (deposit)
# Output: Figure1.pdf, Figure1.png, and printed panel values.
# Usage:  Rscript scripts/figures/make_Fig1.R

suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2)
  library(UpSetR); library(cowplot); library(grid); library(png) })
ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

sp   <- c("col0","ost0","aspen","birch","spruce","pine")
labs <- c(col0="Col-0", ost0="Ost-0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")
key  <- c(col0="col0", ost0="ost0", aspen="aspen", birch="birch", spruce="spruce", pine="pine")

# gene -> orthogroup maps (objects OG_Col0, OG_Ost0, OG_aspen, OG_birch, OG_spruce, OG_pine)
for (f in c("Col0","Ost0","aspen","birch","spruce","pine"))
  load(sprintf("data/ComPlEx/gene2OG/gene2OG_%s.rda", f))
OGm <- list(col0=OG_Col0, ost0=OG_Ost0, aspen=OG_aspen, birch=OG_birch, spruce=OG_spruce, pine=OG_pine)

# ── panel a: expressed-orthogroup sets per condition + UpSet ──────────────────
ctl <- list(); trt <- list()
for (s in sp) {
  e <- suppressMessages(as.data.frame(readr::read_delim(
        sprintf("data/expression/med_expression/%s_med.txt.gz", s), show_col_types = FALSE)))
  cc <- grep("Control$", colnames(e), value = TRUE)
  tc <- grep("(6h|24h|3d|10d)$", colnames(e), value = TRUE)
  ec <- e$Genes[ rowSums(as.matrix(e[, cc, drop = FALSE]) > 0) == length(cc) ]
  et <- e$Genes[ rowSums(as.matrix(e[, tc, drop = FALSE]) > 0) > 0 ]
  og <- OGm[[s]]; gc <- key[s]
  ctl[[s]] <- unique(og$OrthoGroup[og[[gc]] %in% ec])
  trt[[s]] <- unique(og$OrthoGroup[og[[gc]] %in% et])
}
sets <- c(setNames(ctl, paste(unname(labs), "control")),
          setNames(trt, paste(unname(labs), "treatment")))
common_OG_expr <- Reduce(intersect, sets)
cat(sprintf("Panel a: all-12 shared expressed orthogroups = %d\n", length(common_OG_expr)))

all_OG <- unique(unlist(sets)); bm <- sapply(sets, function(x) all_OG %in% x)
OG_binary <- data.frame(lapply(as.data.frame(bm), as.numeric), row.names = all_OG, check.names = FALSE)

# ── panel b: expressed genes in the shared core, per species ──────────────────
bg <- list(col0 = ld("data/annotation/bg_col0_ATC.rda"), ost0 = ld("data/annotation/bg_ost0_ATO.rda"),
           aspen = ld("data/annotation/bg_aspen_PT.rda"), birch = ld("data/annotation/bg_birch_BP.rda"),
           spruce = ld("data/annotation/bg_spruce_PA.rda"), pine = ld("data/annotation/bg_pine_PS.rda"))
tile <- function(df, title, pal) ggplot(df, aes(1, Percent, fill = Group)) +
  geom_col(width = 0.5, colour = "white") +
  geom_text(aes(label = paste0(Group, "\n", Count, " (", round(Percent, 1), "%)")),
            position = position_stack(vjust = 0.5), size = 3.4) +
  scale_fill_manual(values = pal) + coord_flip() + labs(title = title) + theme_void() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5), legend.position = "none")
pal6 <- rev(c("#7CC347","#95EC96","#DAA918","#F3D476","#A770C5","#DEADF9"))
Bcount <- sapply(sp, function(s){ og <- OGm[[s]]; gc <- key[s]
  length(unique(og[[gc]][ og[[gc]] %in% bg[[s]] & og$OrthoGroup %in% common_OG_expr ])) })
Bdf <- data.frame(Group = factor(unname(labs), levels = rev(unname(labs))),
                  Count = as.integer(Bcount), Percent = 100 * Bcount / sapply(bg, length))
pB <- tile(Bdf, "Composition of expressed genes in common orthogroups", pal6)

# ── panels c/d: DEG categories ────────────────────────────────────────────────
DEGl  <- lapply(sp, function(s) ld(sprintf("data/DEGs/DEGs_%s.RData", s))); names(DEGl) <- sp
degOG <- lapply(sp, function(s){ og <- OGm[[s]]; unique(og$OrthoGroup[og[[key[s]]] %in% DEGl[[s]]]) }); names(degOG) <- sp
allOGd <- unique(unlist(degOG)); nsp <- sapply(allOGd, function(o) sum(sapply(degOG, function(x) o %in% x))); names(nsp) <- allOGd
core_OG_DEG <- names(nsp[nsp == 6])
cat(sprintf("Panels c/d: orthogroups with a DEG in all six species = %d\n", length(core_OG_DEG)))

catlev <- c("DEG across all","DEG in some","Species-specific DEG","Singletons")
Crows <- list(); Dcount <- integer(0)
for (s in sp) { og <- OGm[[s]]; D <- DEGl[[s]]; tot <- length(D)
  ogof <- og$OrthoGroup[match(D, og[[key[s]]])]; nc <- nsp[ogof]
  cnts <- c("DEG across all" = sum(!is.na(nc) & nc == 6), "DEG in some" = sum(!is.na(nc) & nc >= 2 & nc <= 5),
            "Species-specific DEG" = sum(!is.na(nc) & nc == 1), "Singletons" = sum(is.na(ogof)))
  Crows[[s]] <- data.frame(species = labs[s], category = names(cnts), Counts = as.integer(cnts),
                           prop = as.numeric(cnts)/tot, row.names = NULL)
  Dcount[s] <- sum(!is.na(ogof) & ogof %in% core_OG_DEG)
}
Cdf <- bind_rows(Crows); Cdf$category <- factor(Cdf$category, levels = catlev)
Cdf$species <- factor(Cdf$species, levels = unname(labs))
pal4 <- rev(c("#A092E4","#80C5BA","#E4A1AC","#CFD676"))
pC <- ggplot(Cdf, aes(species, Counts, fill = category)) + geom_col(colour = "black", width = 0.7) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = pal4) + theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Number of DEGs", fill = "Category") +
  theme(panel.grid = element_blank(), legend.position = "top")
Ddf <- data.frame(Group = factor(unname(labs), levels = rev(unname(labs))),
                  Count = as.integer(Dcount[sp]), Percent = 100 * Dcount[sp] / sapply(DEGl, length))
pD <- tile(Ddf, "Composition of DEGs in common orthogroups", pal6)

cat("\nPanel b (expressed genes in common orthogroups):\n")
for (s in sp) cat(sprintf("  %-6s %d (%.1f%%)\n", labs[s], Bdf$Count[Bdf$Group==labs[s]], Bdf$Percent[Bdf$Group==labs[s]]))
cat("Panel d (DEGs in common orthogroups):\n")
for (s in sp) cat(sprintf("  %-6s %d (%.1f%%)\n", labs[s], Ddf$Count[Ddf$Group==labs[s]], Ddf$Percent[Ddf$Group==labs[s]]))

# ── render ────────────────────────────────────────────────────────────────────
png("panelA_fig1.png", width = 2600, height = 1150, res = 210)
upset(OG_binary, sets = rev(colnames(OG_binary)), keep.order = TRUE, nintersects = 15,
      order.by = "freq", decreasing = TRUE, set_size.show = TRUE,
      main.bar.color = "#636363", matrix.color = "#828282", sets.bar.color = "#828282",
      sets.x.label = "Total orthogroups per selection", mainbar.y.label = "Orthogroup overlap of expressed genes",
      set_size.scale_max = 18000, text.scale = c(1.5,1.5,0,1.5,1.5,1.6))
dev.off()
pA <- ggdraw() + draw_image("panelA_fig1.png")
top <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1.9, 1), labels = c("a","b"))
bot <- plot_grid(pC, pD, ncol = 2, rel_widths = c(1.4, 1), labels = c("c","d"))
fig <- plot_grid(top, bot, ncol = 1, rel_heights = c(1.25, 1))
ggsave("Figure1.pdf", fig, width = 15, height = 11)
ggsave("Figure1.png", fig, width = 15, height = 11, dpi = 110)
cat("\nWrote Figure1.pdf / Figure1.png\n")
