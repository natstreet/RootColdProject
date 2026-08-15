#!/usr/bin/env Rscript
# make_FigS3.R — Figure S3: for each species, the number of genes classed as co-expressologs and
# as differentially expressed (DE) co-expressologs, split by whether the partner gene in the other
# species is itself differentially expressed. Co-expressologs are taken at MaxpVal < 0.1.
#
# Run from a directory containing the unpacked figshare deposit `data/`.
# Inputs: data/ComPlEx/co_expressologs.RData, data/DEGs/DEGs_<species>.RData
# Output: Figure_S3.pdf, Figure_S3.png
# Usage:  Rscript scripts/figures/make_FigS3.R

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
P <- 0.1
specs <- c("col0","ost0","aspen","birch","spruce","pine")
labmap <- c(col0="Col-0", ost0="Ost-0", aspen="Aspen", birch="Birch", spruce="Spruce", pine="Pine")

ld <- function(f){ e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
co <- ld("data/ComPlEx/co_expressologs.RData")
DEG <- lapply(specs, function(s) as.character(ld(sprintf("data/DEGs/DEGs_%s.RData", s)))); names(DEG) <- specs

# each co-expressolog edge contributes two (gene, species, partner) rows
L <- bind_rows(
  transmute(co, gene = GeneSpecies1, sp = Species1, pg = GeneSpecies2, psp = Species2, MaxpVal),
  transmute(co, gene = GeneSpecies2, sp = Species2, pg = GeneSpecies1, psp = Species1, MaxpVal)) %>%
  filter(MaxpVal < P)

cat_counts <- function(S){
  d <- L %>% filter(sp == S)
  d$isDE   <- d$gene %in% DEG[[S]]
  d$partDE <- mapply(function(pg, psp) pg %in% DEG[[psp]], d$pg, d$psp)
  g <- d %>% group_by(gene) %>% summarise(isDE = any(isDE), partDE = any(partDE), .groups = "drop")
  c(coexpr = sum(!g$isDE & !g$partDE), coexpr_DEpart = sum(!g$isDE & g$partDE),
    DEcoexpr = sum(g$isDE & !g$partDE), DEcoexpr_DEpart = sum(g$isDE & g$partDE))
}
m <- as.data.frame(t(sapply(specs, cat_counts))); m$species <- labmap[rownames(m)]

catlev <- c("DEcoexpr_DEpart","DEcoexpr","coexpr_DEpart","coexpr")
catlab <- c(DEcoexpr_DEpart = "DE co-expressologs, partner DE in another species",
            DEcoexpr        = "DE co-expressologs",
            coexpr_DEpart   = "co-expressologs, partner DE in another species",
            coexpr          = "co-expressologs")
long <- m %>% pivot_longer(all_of(catlev), names_to = "cat", values_to = "n")
long$species <- factor(long$species, levels = unname(labmap))
long$cat <- factor(long$cat, levels = catlev, labels = catlab[catlev])

pal <- setNames(c("#1B5E9C","#4A90C2","#C46E00","#E0B080"), catlab[catlev])
p <- ggplot(long, aes(species, n, fill = cat)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = NULL, y = "Number of genes") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), legend.position = "top",
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 2))
ggsave("Figure_S3.pdf", p, width = 9, height = 6)
ggsave("Figure_S3.png", p, width = 9, height = 6, dpi = 150)
print(m[, c("species", catlev)], row.names = FALSE)
