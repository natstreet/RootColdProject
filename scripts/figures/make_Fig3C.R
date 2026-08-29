#!/usr/bin/env Rscript
# make_Fig3C.R — Figure 3C: conservation of each supercluster's genes across the three
# phylogenetic groups (Arabidopsis / Broadleaf / Conifer), for DEGs and TF-DEGs.
#
# Method: for each of the four superclusters carried into Figure 3, map every member gene to
# its ortholog group (OrthoFinder), and record which of the three groups contribute a gene to
# that supercluster. The panel shows the distribution of the supercluster's unique orthogroups
# over the seven categories (three single-group, three pairwise, all-three) as 100% stacked bars.
#
# The main-text labels SC_I-SC_IV are a renumbering of four of the twelve superclusters in
# Fig. S1; the correspondence is read from SC_mapping.csv (SC_I=SC6, SC_II=SC2, SC_III=SC3,
# SC_IV=SC1).
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/superclusters/SC_mapping.csv
#   data/superclusters/DEGs_TFs/SC<n>/SC<n>_gene_list.csv
#   data/ComPlEx/Orthogroups_20240613.tsv   (or data/annotation/…)
#   data/TranscriptionFactors/TF_<species>.RData   (for the TF-DEGs)
# Outputs:
#   figures/Figure3C.png ; figures/Fig3C_values.csv
#
# Usage: Rscript scripts/figures/make_Fig3C.R

suppressMessages({ library(ggplot2) })

ogfile <- Sys.glob("data/*/Orthogroups_20240613.tsv")[1]
og <- read.delim(ogfile, stringsAsFactors = FALSE, check.names = FALSE); cols <- colnames(og)
Ar <- "Arabidopsis_thaliana"; Bi <- "Betula_pendula"
Po <- grep("Populus_trem", cols, value = TRUE)[1]
Pa <- grep("Picea_abies", cols, value = TRUE)[1]; Ps <- grep("Pinus_sylv", cols, value = TRUE)[1]
splt <- function(x) unlist(strsplit(x, "[, ]+"))
g2og <- new.env(hash = TRUE)
for (i in seq_len(nrow(og)))
  for (g in c(splt(og[i, Ar]), splt(og[i, Bi]), splt(og[i, Po]), splt(og[i, Pa]), splt(og[i, Ps])))
    if (nzchar(g)) assign(g, i, envir = g2og)
clade <- function(g) if (grepl("^AT", g)) "A" else if (grepl("^Bpev|^Potr", g)) "D" else "C"

getobj <- function(f) { e <- new.env(); load(f, envir = e); x <- get(ls(e)[1], envir = e); if (is.data.frame(x)) x[[1]] else x }
TF <- unique(sub("\\.\\d+$", "", unlist(lapply(Sys.glob("data/TranscriptionFactors/TF_*.RData"), getobj))))

mapping <- read.csv("data/superclusters/SC_mapping.csv", stringsAsFactors = FALSE)  # Main_text, Original_supercluster
CAT <- c("A","C","D","AC","AD","CD","ACD")
NAMES <- c(A="Arabidopsis only", C="Conifer only", D="Broadleaf only",
           AC="Arabidopsis & Conifer", AD="Arabidopsis & Broadleaf",
           CD="Conifer & Broadleaf", ACD="All three groups")
share_of <- function(orig, tfonly) {
  f <- Sys.glob(sprintf("data/superclusters/DEGs_TFs/%s/%s_gene_list.csv", orig, orig))
  sc <- sub("\\.\\d+$", "", read.csv(f[1], stringsAsFactors = FALSE)[, 2])
  if (tfonly) sc <- sc[sc %in% TF]
  ogi <- sapply(sc, function(g) { i <- mget(g, envir = g2og, ifnotfound = NA)[[1]]; if (is.null(i) || is.na(i)) NA else i })
  cl <- sapply(sc, clade); ok <- !is.na(ogi)
  ogset <- tapply(cl[ok], ogi[ok], function(v) paste(sort(unique(v)), collapse = ""))
  t <- table(factor(ogset, levels = CAT)); 100 * t / sum(t)
}

res <- list()
for (k in seq_len(nrow(mapping))) for (ty in c("DEGs","TF-DEGs")) {
  v <- share_of(mapping$Original_supercluster[k], ty == "TF-DEGs")
  res[[length(res)+1]] <- data.frame(SC = mapping$Main_text[k], type = ty,
                                     category = factor(NAMES[CAT], levels = NAMES[CAT]), pct = as.numeric(v))
}
D <- do.call(rbind, res)
dir.create("figures", showWarnings = FALSE)
write.csv(reshape(transform(D, category = names(NAMES)[match(category, NAMES)]),
                  idvar = c("SC","type"), timevar = "category", direction = "wide"),
          "figures/Fig3C_values.csv", row.names = FALSE)

# 100% stacked bars, categories ordered lineage-specific -> shared by two -> all three
D$category <- factor(D$category, levels = NAMES[c("A","D","C","AD","AC","CD","ACD")])
D$SC   <- factor(D$SC, levels = rev(mapping$Main_text))
D$type <- factor(D$type, levels = c("DEGs","TF-DEGs"))
pal <- c("Arabidopsis only"="#2F9E6E", "Broadleaf only"="#D8892B", "Conifer only"="#7B6BA8",
         "Arabidopsis & Broadleaf"="#6FA88C", "Arabidopsis & Conifer"="#5E8CA8",
         "Conifer & Broadleaf"="#A67BA0", "All three groups"="#4D4D4D")
p <- ggplot(D, aes(pct, SC, fill = category)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct >= 6, round(pct), ""), group = category),
            position = position_stack(vjust = 0.5), colour = "white", size = 2.7) +
  facet_wrap(~ type, ncol = 1, strip.position = "top") +
  scale_fill_manual(values = pal, name = NULL, guide = guide_legend(ncol = 1)) +
  scale_x_continuous(expand = expansion(0), breaks = c(0, 25, 50, 75, 100)) +
  labs(x = "Orthogroups (%)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_text(face = "bold", colour = "black"), axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0), strip.background = element_blank(),
        legend.position = "right")
ggsave("figures/Figure3C.png", p, width = 8, height = 6, dpi = 120)
cat("Wrote figures/Figure3C.png and figures/Fig3C_values.csv\n")
print(reshape(D, idvar = c("SC","type"), timevar = "category", direction = "wide"), row.names = FALSE)
