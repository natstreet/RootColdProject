#!/usr/bin/env Rscript
# make_Fig6.R — Figure 6: expression profiles across species of transcription factors from the
# conserved co-expression core. One panel per transcription factor (labelled by its Arabidopsis
# symbol); within a panel, the VST expression profile of that orthogroup's member in each species
# over the cold time course (mean line, +/-1 SD ribbon, individual replicates as points). A
# coloured star marks each species in which that member is differentially expressed (padj < 0.05).
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs:
#   data/annotation/gene_aliases_20140331.txt              Arabidopsis symbol -> gene ID
#   data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS  orthogroup -> per-species member
#   data/expression/<species>_root_cold_stress_expression.txt  per-sample VST
#   data/DEGs/DEGs_<species>.RData                         per-species DEG sets (for the stars)
# Output: Figure6.pdf, Figure6.png
#
# Usage: Rscript scripts/figures/make_Fig6.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

symbols <- c("ANAC018","bZIP61","HSFB2A","KAN3","LBD11","LBD13","MYB87","WRKY48")  # panel order
species <- c(col0="Col-0", ost0="Ost-0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")
pre     <- c(col0="Col", ost0="Ost", aspen="Pt", birch="Bp", spruce="Pa", pine="Ps")
cond    <- c("01"="Control","02"="6h","03"="24h","04"="3d","05"="10d")
pal     <- c("Col-0"="#2ca02c","Ost-0"="#ff7f0e","Aspen"="#1f3b8c","Birch"="#4aa3df",
             "Spruce"="#d62728","Pine"="#e377c2")
ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

# symbol -> Arabidopsis gene ID
al  <- read.delim("data/annotation/gene_aliases_20140331.txt", stringsAsFactors = FALSE)
sym2id <- function(s){ i <- which(toupper(al$symbol) == toupper(s)); if (!length(i)) NA else al$locus_name[i[1]] }
id2sym <- function(g){ i <- which(al$locus_name == g); if (!length(i)) g else al$symbol[i[1]] }
seed <- setNames(vapply(symbols, sym2id, character(1)), symbols)

# clique table -> orthogroup and its per-species members
cl <- readRDS("data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS")
og_of <- function(gid){ h <- cl[cl$GeneSpecies1 == gid | cl$GeneSpecies2 == gid, "OrthoGroup", drop = TRUE]; if (!length(h)) NA else h[1] }
members <- function(og){                                   # species -> gene id(s)
  r <- cl[cl$OrthoGroup == og, ]
  m <- tapply(c(r$GeneSpecies1, r$GeneSpecies2), c(r$Species1, r$Species2), function(x) unique(x))
  lapply(m, as.character)
}

# per-species per-sample VST for one gene -> long table of (Species, Condition, rep, VST)
DEG <- lapply(names(species), function(sp) as.character(ld(sprintf("data/DEGs/DEGs_%s.RData", sp))))
names(DEG) <- names(species)
EXP <- lapply(names(species), function(sp)
  read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp), check.names = FALSE, row.names = 1))
names(EXP) <- names(species)

rows <- list(); stars <- list()
for (s in symbols) {
  og <- og_of(seed[[s]]); if (is.na(og)) next
  mem <- members(og)
  lab <- s                                                # panel title = the Arabidopsis symbol
  for (sp in names(species)) {
    g <- intersect(mem[[sp]], rownames(EXP[[sp]])); if (!length(g)) next
    v <- colMeans(EXP[[sp]][g, , drop = FALSE], na.rm = TRUE)                          # if >1 paralog, average
    cc <- cond[sub(paste0("^", pre[sp], "\\d+\\."), "", names(v))]
    df <- data.frame(TF = lab, Species = species[[sp]], Condition = cc, VST = as.numeric(v))
    rows[[paste(s, sp)]] <- df
    if (any(g %in% DEG[[sp]])) stars[[paste(s, sp)]] <- data.frame(TF = lab, Species = species[[sp]])
  }
}
d <- bind_rows(rows)
d$Condition <- factor(d$Condition, levels = unname(cond))
d$Species   <- factor(d$Species,   levels = unname(species))
d$TF        <- factor(d$TF, levels = unique(d$TF))
summ <- d %>% group_by(TF, Species, Condition) %>%
  summarise(m = mean(VST), sd = sd(VST), .groups = "drop")

# place the DEG stars just below the top of each panel, stacked per species
star <- bind_rows(stars)
if (nrow(star)) {
  star$TF <- factor(star$TF, levels = levels(d$TF)); star$Species <- factor(star$Species, levels = levels(d$Species))
  top <- d %>% group_by(TF) %>% summarise(ymax = max(VST), .groups = "drop")
  star <- star %>% left_join(top, by = "TF") %>% group_by(TF) %>%
    mutate(y = ymax - (row_number() - 1) * 0.06 * ymax) %>% ungroup()
}

p <- ggplot(summ, aes(Condition, m, colour = Species, group = Species)) +
  geom_ribbon(aes(ymin = m - sd, ymax = m + sd, fill = Species), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(data = d, aes(Condition, VST, colour = Species), size = 0.5, alpha = 0.5) +
  { if (nrow(star)) geom_text(data = star, aes(x = length(cond) + 0.4, y = y, colour = Species),
        label = "*", size = 5, show.legend = FALSE) } +
  facet_wrap(~TF, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal, guide = "none") +
  labs(x = "Condition", y = "Expression (VST)", colour = "Species",
       title = "Species expression profiles of transcription factors",
       caption = "* differentially expressed in that species (DEG, padj < 0.05)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figure6.pdf", p, width = 13, height = 6)
ggsave("Figure6.png", p, width = 13, height = 6, dpi = 150)
cat("TFs plotted:", nlevels(d$TF), "| orthogroups:", paste(levels(d$TF), collapse=", "), "\n")
