#!/usr/bin/env Rscript
# make_Fig2.R — Figure 2: per-timepoint numbers of up- and down-regulated DEGs per species.
#
# Reproduces Figure 2 of Aro et al. from the deposited differential-expression results.
#
# DEG definition: adjusted P <= 0.05 AND |log2 fold change| >= 1 (a two-fold change). The
# `..._p0.05_lfc0.txt` per-timepoint lists used here are the DESeq2 Wald-test contrasts of each
# cold timepoint against the pooled control; every gene in them satisfies |log2FC| >= 1.
# Direction (up/down) is the sign of the log2 fold change.
#
# Inputs  (unpack data.tar.gz from the figshare deposit into ./data):
#   data/DEGs/per_timepoint_DEG_lists/<species>/DEGs_<sp>_C_vs_<tp>_p0.05_lfc0.txt
#   data/DEGs/log2FoldChange_<code>.RData   (per-timepoint log2FC matrices; for direction)
# Output:
#   figures/Figure2.png
#
# Usage:  Rscript scripts/figures/make_Fig2.R

suppressMessages(library(ggplot2))

tps  <- c("6h", "24h", "3d", "10d")
# species -> (DEG-list filename code, log2FoldChange matrix code)
sp <- list(
  "Col-0 roots"  = c(list = "col0", lfc = "col0"),
  "Ost-0 roots"  = c(list = "ost0", lfc = "ost0"),
  "Aspen roots"  = c(list = "pt",   lfc = "pt"),
  "Birch roots"  = c(list = "bp",   lfc = "bp"),
  "Spruce roots" = c(list = "pa",   lfc = "pa"),
  "Pine roots"   = c(list = "ps",   lfc = "ps"))
dirs <- c("Col-0 roots" = "col0", "Ost-0 roots" = "ost0", "Aspen roots" = "aspen",
          "Birch roots" = "birch", "Spruce roots" = "spruce", "Pine roots" = "pine")

getobj <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

rows <- list()
for (lab in names(sp)) {
  code <- sp[[lab]]["list"]; lfc <- getobj(sprintf("data/DEGs/log2FoldChange_%s.RData", sp[[lab]]["lfc"]))
  rn <- rownames(lfc)
  for (tp in tps) {
    f <- Sys.glob(sprintf("data/DEGs/per_timepoint_DEG_lists/%s/DEGs_%s_C_vs_%s_p0.05_lfc0.txt",
                          dirs[lab], code, tp))
    genes <- readLines(f[1]); genes <- genes[nzchar(genes)]
    v <- setNames(lfc[[paste0("log2FoldChange_", tp)]], rn)[genes]; v <- v[!is.na(v)]
    rows[[length(rows) + 1]] <- data.frame(species = lab, tp = tp,
      Direction = c("Up", "Down"), Count = c(sum(v > 0), sum(v < 0)))
  }
}
d <- do.call(rbind, rows)
d$species   <- factor(d$species, levels = names(sp))
d$tp        <- factor(d$tp, levels = tps)
d$Direction <- factor(d$Direction, levels = c("Up", "Down"))
# tag panels (A)-(F)
levels(d$species) <- paste0("(", LETTERS[seq_along(levels(d$species))], ") ", levels(d$species))

p <- ggplot(d, aes(tp, Count, fill = Direction)) +
  geom_col(width = 0.8) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Up = "#C03000", Down = "#00B0F0")) +
  labs(x = "Timepoint", y = "Number of DEGs") +
  theme_bw(base_size = 16) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 15, hjust = 0),
        axis.text = element_text(colour = "black"))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure2.pdf", p, width = 12.54, height = 12.54)
ggsave("figures/Figure2.png", p, width = 12.54, height = 12.54, dpi = 110)
cat("Wrote figures/Figure2.pdf and figures/Figure2.png\n")
