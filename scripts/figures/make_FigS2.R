#!/usr/bin/env Rscript
# make_FigS2.R — Figure S2: transcription-factor family composition of each super cluster,
# resolved per species. Each panel is one super cluster (SC_I–SC_IV); within a panel, one
# stacked bar per species/ecotype gives the percentage of that super cluster's transcription
# factors falling in each of the five most represented families (ERF, MYB, NAC, bHLH, WRKY),
# with all remaining families pooled as "Others".
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs:
#   data/superclusters/membership/<species>/cluster_<n>_members.csv
#       per-species gene membership of each super cluster (original numbering); the
#       transcription factors among them are selected via the TF family table below
#   data/TranscriptionFactors/TF_<species>.RData   gene -> TF family (Ost-0 uses the Col-0 set)
# Output: Figure_S2.pdf, Figure_S2.tiff
#
# Usage: Rscript scripts/figures/make_FigS2.R

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })

# main-text super cluster label -> original super cluster number
SC <- c(SC_I = 6, SC_II = 2, SC_III = 3, SC_IV = 1)
species <- c(col0 = "Col-0", ost0 = "Ost-0", birch = "Birch",
             aspen = "Aspen", pine = "Pine", spruce = "Spruce")
majors  <- c("ERF", "MYB", "NAC", "bHLH", "WRKY")

ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
famtab <- function(sp){
  x <- ld(sprintf("data/TranscriptionFactors/TF_%s.RData", if (sp == "ost0") "col0" else sp))
  if (all(c("Gene_ID","Family") %in% colnames(x)))
    setNames(as.character(x$Family), as.character(x$Gene_ID))
  else setNames(as.character(x[[2]]), as.character(x[[1]]))
}
fam <- lapply(names(species), famtab); names(fam) <- names(species)

rows <- list()
for (lab in names(SC)) for (sp in names(species)) {
  f <- sprintf("data/superclusters/membership/%s/cluster_%d_members.csv", sp, SC[[lab]])
  g <- read.csv(f, stringsAsFactors = FALSE)$gene_id
  fv <- fam[[sp]][as.character(g)]; fv <- fv[!is.na(fv)]
  cls <- ifelse(fv %in% majors, fv, "Others")
  tb  <- table(factor(cls, levels = c(majors, "Others")))
  rows[[paste(lab, sp)]] <- data.frame(
    SC = factor(lab, levels = names(SC)), Species = species[[sp]],
    Family = names(tb), n = as.integer(tb), total = sum(tb))
}
df <- bind_rows(rows) %>% mutate(pct = ifelse(total > 0, 100 * n / total, 0))
df$Species <- factor(df$Species, levels = unname(species))
df$Family  <- factor(df$Family, levels = c("Others", "WRKY", "bHLH", "NAC", "MYB", "ERF"))

n_lab <- unique(df[, c("SC", "Species", "total")])          # transcription factors per species/SC
pal <- c(ERF = "#EEC200", MYB = "#3B4CB8", NAC = "#9C3A2E",
         bHLH = "#111111", WRKY = "#4A235A", Others = "#8C8C8C")
p <- ggplot(df, aes(Species, pct, fill = Family)) +
  geom_col(width = 0.8) +
  geom_text(data = n_lab, aes(Species, 100, label = total), inherit.aes = FALSE,
            vjust = -0.4, size = 2.6, colour = "grey30") +
  facet_wrap(~SC, ncol = 2) +
  scale_fill_manual(values = pal, breaks = c(majors, "Others"), name = "TF family") +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0.09))) +
  labs(x = NULL, y = "Percentage (%)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey92"))

ggsave("Figure_S2.pdf",  p, width = 9.5, height = 6)
ggsave("Figure_S2.tiff", p, width = 9.5, height = 6, dpi = 300, compression = "lzw")
