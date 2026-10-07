#!/usr/bin/env Rscript
# make_Fig6.R — Figure 6: expression profiles across species of transcription factors from the
# conserved co-expression core. One panel per transcription factor (labelled by its Arabidopsis
# symbol); within a panel, the VST expression profile of that orthogroup's member in each species
# over the cold time course (mean line, +/-1 SD ribbon, individual replicates as points). A
# coloured star marks each species in which that member is differentially expressed (adjusted P <= 0.05
# and |log2FC| >= 1 at any timepoint).
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs:
#   scripts/analysis/tf_cliques.R (the transcription-factor cliques: >= 5 of 6 members predicted TFs, and their
#   labels), which reads data/TranscriptionFactors/ and data/annotation/gene_aliases_20140331.txt
#   data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS orthogroup -> per-species clique member
#   data/expression/<species>_root_cold_stress_expression.txt  per-sample VST
#   data/DEGs/DEGs_<species>.RData                         per-species DEG sets (for the stars)
#   results/orthogroup_size_flag.tsv (scripts/validation/orthogroup_size_flag.R): each panel's label carries its
#   orthogroup's number of expressed genes (NRS, 2026-10-05)
#   results/orthology_class.tsv (scripts/analysis/orthology_class.R): and its orthology class (D043)
# Output: Figure6.pdf, Figure6.png
#
# Usage: Rscript scripts/figures/make_Fig6.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

source("scripts/analysis/tf_cliques.R")
TFC <- tf_cliques()                       # the transcription-factor cliques of the core, labelled, in panel order
species <- c(col0="Col-0", ost0="Ost-0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")
pre     <- c(col0="Col", ost0="Ost", aspen="Pt", birch="Bp", spruce="Pa", pine="Ps")
cond    <- c("01"="Control","02"="6h","03"="24h","04"="3d","05"="10d")
source("scripts/figures/genotype_colours.R")
pal     <- setNames(GENOTYPE_COLOURS[names(species)], species)   # the shared genotype palette, by display name
ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

# the selected conserved clique: one member gene per species for each orthogroup
cgsel <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
members <- function(og){                                   # species -> that species' clique member
  r <- cgsel[cgsel$OrthoGroup == og, ]
  if (!nrow(r)) return(list())
  setNames(lapply(names(species), function(sp) as.character(r[[sp]])), names(species))
}

# per-species per-sample VST for one gene -> long table of (Species, Condition, rep, VST)
DEG <- lapply(names(species), function(sp) as.character(ld(sprintf("data/DEGs/DEGs_%s.RData", sp))))
names(DEG) <- names(species)
EXP <- lapply(names(species), function(sp)
  read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", sp), check.names = FALSE, row.names = 1))
names(EXP) <- names(species)

SIZE <- read.delim("results/orthogroup_size_flag.tsv")
if (!all(TFC$OrthoGroup %in% SIZE$OrthoGroup)) stop("a transcription-factor clique has no size flag")
OCL <- read.delim("results/orthology_class.tsv")
if (!all(TFC$OrthoGroup %in% OCL$OrthoGroup)) stop("a transcription-factor clique has no orthology class")
rows <- list(); stars <- list()
for (k in seq_len(nrow(TFC))) {
  og <- TFC$OrthoGroup[k]; s <- TFC$tf[k]                  # TF symbol (dual where Ost-0 member differs)
  lab <- sprintf("%s (%d genes, %s)", s, SIZE$expressed_genes[SIZE$OrthoGroup == og],
                 OCL$orthology_class[OCL$OrthoGroup == og])                           # panel label
  mem <- members(og)
  for (sp in names(species)) {
    g <- intersect(mem[[sp]], rownames(EXP[[sp]])); if (!length(g)) next
    v <- colMeans(EXP[[sp]][g, , drop = FALSE], na.rm = TRUE)                          # the clique member
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
       caption = expression("* differentially expressed in that species (DEG: adjusted" ~ italic(P) <= 0.05 ~
                            "and |" * log[2] * "FC|" >= 1 ~ "at any timepoint)")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figure6.pdf", p, width = 13, height = 6, bg = "white")
ggsave("Figure6.png", p, width = 13, height = 6, dpi = 300, bg = "white")
if (nlevels(d$TF) != nrow(TFC)) stop("plotted ", nlevels(d$TF), " of ", nrow(TFC), " transcription-factor cliques")
cat("TFs plotted:", nlevels(d$TF), "| orthogroups:", paste(levels(d$TF), collapse=", "), "\n")
