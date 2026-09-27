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
#   data/DEGs/og_summary.RData                               per-orthogroup species counts (deposit)
# Output: Figure1.pdf, Figure1.png, and Figure1_panel_values.tsv (panels a-d: count, denominator, percent).
# Usage:  Rscript scripts/figures/make_Fig1.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2)
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
# one equal-width tile per genotype, labelled with its count and percentage
tile <- function(df, title, pal) {
  df$x <- match(as.character(df$Group), unname(labs[sp]))
  ggplot(df, aes(x, 1, fill = Group)) +
    geom_tile(width = 0.98, height = 1) +
    geom_text(aes(label = paste0(Group, "\n", Count, "\n(", round(Percent, 1), "%)")), size = 3.6, lineheight = 0.95) +
    scale_fill_manual(values = pal) + scale_x_continuous(expand = c(0, 0)) + labs(title = title) + theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 6)),
          legend.position = "none")
}
source("scripts/figures/genotype_colours.R")
pal6 <- setNames(genotype_tint(GENOTYPE_COLOURS[sp]), labs[sp])   # the shared genotype palette, lightened for black labels
Bcount <- sapply(sp, function(s){ og <- OGm[[s]]; gc <- key[s]
  length(unique(og[[gc]][ og[[gc]] %in% bg[[s]] & og$OrthoGroup %in% common_OG_expr ])) })
Bdf <- data.frame(Group = factor(unname(labs), levels = rev(unname(labs))),
                  Count = as.integer(Bcount), Percent = 100 * Bcount / sapply(bg, length))
pB <- tile(Bdf, "Composition of expressed genes in common orthogroups", pal6)

# ── panels c/d: DEG categories (deposited orthogroup table + og_summary) ──────
# Each DEG is mapped to its orthogroup through the per-species gene -> orthogroup maps (OGm) and
# classified by how many species carry a gene / a DEG in that orthogroup (deposited og_summary).
DEGl  <- lapply(sp, function(s) ld(sprintf("data/DEGs/DEGs_%s.RData", s))); names(DEGl) <- sp
orthos <- bind_rows(lapply(sp, function(s)
  data.frame(OrthoGroup = OGm[[s]]$OrthoGroup, GeneID = OGm[[s]][[key[s]]], species = s,
             stringsAsFactors = FALSE)))
orthos$DEG <- mapply(function(g, s) g %in% DEGl[[s]], orthos$GeneID, orthos$species)
ogsum <- ld("data/DEGs/og_summary.RData")   # per-orthogroup n_species_with_genes / n_species_with_DEGs
DEGclasses <- orthos %>% filter(DEG) %>% left_join(ogsum, by = "OrthoGroup") %>%
  mutate(category = case_when(
    n_species_with_genes == 1 & n_species_with_DEGs == 1 ~ "Singletons",
    n_species_with_DEGs == 6                             ~ "DEG across all",
    n_species_with_DEGs >  1 & n_species_with_DEGs != 6  ~ "DEG in some",
    n_species_with_DEGs == 1 & n_species_with_genes >  1 ~ "Species-specific DEG",
    TRUE ~ NA_character_))
cat(sprintf("Panels c/d: orthogroups with a DEG in all six species = %d\n", sum(ogsum$n_species_with_DEGs == 6)))

catlev <- c("DEG across all","DEG in some","Species-specific DEG","Singletons")
Crows <- list(); Dcount <- integer(0)
for (s in sp) {
  sub  <- DEGclasses[DEGclasses$species == s & !is.na(DEGclasses$category), ]
  cnts <- setNames(sapply(catlev, function(c) sum(sub$category == c)), catlev)
  Crows[[s]] <- data.frame(species = labs[s], category = catlev, Counts = as.integer(cnts),
                           prop = as.numeric(cnts)/sum(cnts), row.names = NULL)
  Dcount[s] <- sum(sub$category == "DEG across all")     # DEGs in orthogroups with a DEG in all six
}
Cdf <- bind_rows(Crows); Cdf$category <- factor(Cdf$category, levels = catlev)
Cdf$species <- factor(Cdf$species, levels = unname(labs))
pal4 <- rev(c("#A092E4","#80C5BA","#E4A1AC","#CFD676"))
pC <- ggplot(Cdf, aes(species, Counts, fill = category)) + geom_col(colour = "black", width = 0.7) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = pal4) + theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Number of DEGs", fill = "Category") +
  theme(panel.grid = element_blank(), legend.position = "right", axis.line = element_line(colour = "black"),
        plot.margin = margin(40, 5, 5, 5))
Ddf <- data.frame(Group = factor(unname(labs), levels = rev(unname(labs))),
                  Count = as.integer(Dcount[sp]), Percent = 100 * Dcount[sp] / sapply(DEGl, length))
pD <- tile(Ddf, "Composition of DEGs in common orthogroups", pal6)

# panel values behind the Results text (rows of Bdf and Ddf are in `sp` order)
vals <- bind_rows(
  data.frame(panel = "a", species = "all", category = "expressed orthogroups shared by all 12 condition sets",
             count = length(common_OG_expr), denominator = length(all_OG), percent = round(100 * length(common_OG_expr) / length(all_OG), 2)),
  data.frame(panel = "b", species = unname(labs[sp]), category = "expressed genes in common orthogroups",
             count = Bdf$Count, denominator = as.integer(sapply(bg, length)[sp]), percent = round(Bdf$Percent, 2)),
  data.frame(panel = "c", species = as.character(Cdf$species), category = as.character(Cdf$category),
             count = Cdf$Counts, denominator = as.integer(ave(Cdf$Counts, Cdf$species, FUN = sum)),
             percent = round(100 * Cdf$prop, 2)),
  data.frame(panel = "d", species = unname(labs[sp]), category = "DEGs in common orthogroups",
             count = Ddf$Count, denominator = as.integer(sapply(DEGl, length)[sp]), percent = round(Ddf$Percent, 2)),
  data.frame(panel = "d/b", species = unname(labs[sp]),
             category = "DEGs in common orthogroups per expressed gene in common orthogroups",
             count = Ddf$Count, denominator = Bdf$Count, percent = round(100 * Ddf$Count / Bdf$Count, 2)),
  data.frame(panel = "c/d", species = "all", category = "orthogroups with a DEG in all six species",
             count = sum(ogsum$n_species_with_DEGs == 6), denominator = nrow(ogsum), percent = NA))
write.table(vals, "Figure1_panel_values.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nPanel b (expressed genes in common orthogroups):\n")
for (s in sp) cat(sprintf("  %-6s %d (%.1f%%)\n", labs[s], Bdf$Count[Bdf$Group==labs[s]], Bdf$Percent[Bdf$Group==labs[s]]))
cat("Panel d (DEGs in common orthogroups):\n")
for (s in sp) cat(sprintf("  %-6s %d (%.1f%%)\n", labs[s], Ddf$Count[Ddf$Group==labs[s]], Ddf$Percent[Ddf$Group==labs[s]]))

# ── render ────────────────────────────────────────────────────────────────────
png("panelA_fig1.png", width = 3900, height = 1750, res = 300)
upset(OG_binary, sets = rev(colnames(OG_binary)), keep.order = TRUE, nintersects = 15,
      order.by = "freq", decreasing = TRUE, set_size.show = TRUE,
      main.bar.color = "#636363", matrix.color = "#828282", sets.bar.color = "#828282",
      sets.x.label = "Total orthogroups per selection", mainbar.y.label = "Orthogroup overlap of expressed genes",
      set_size.scale_max = 18000, text.scale = c(1.5,1.5,0,1.5,1.5,1.6))
dev.off()
# layout as in the article: panels b and d sit in the empty upper right of panels a and c (13 x 12 in)
fig <- ggdraw() +
  draw_image("panelA_fig1.png", x = 0, y = 0.52, width = 1, height = 0.48) +
  draw_plot(pB, x = 0.47, y = 0.86, width = 0.53, height = 0.13) +
  draw_plot(pC, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot(pD, x = 0.47, y = 0.36, width = 0.53, height = 0.13) +
  draw_plot_label(c("a", "b", "c", "d"), x = c(0, 0.45, 0, 0.45), y = c(1, 1, 0.5, 0.5), size = 16)
ggsave("Figure1.pdf", fig, width = 13, height = 12, bg = "white")
ggsave("Figure1.png", fig, width = 13, height = 12, dpi = 300, bg = "white")
cat("\nWrote Figure1.pdf / Figure1.png\n")
