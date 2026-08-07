#!/usr/bin/env Rscript
# make_Fig3C.R — Figure 3C: conservation of each supercluster's genes across the three
# phylogenetic groups (Arabidopsis / Deciduous / Coniferous), for DEGs and TF-DEGs.
#
# Method: for each of the four superclusters carried into Figure 3, map every member gene to
# its ortholog group (OrthoFinder), and record which of the three groups contribute a gene to
# that supercluster. The pie is the distribution of the supercluster's unique orthogroups over
# the seven categories (three single-group, three pairwise, all-three).
#
# The main-text labels SC_I-SC_IV are a renumbering of four of the twelve superclusters in
# Fig. S1; the correspondence is read from SC_mapping.csv (SC_I=SC6, SC_II=SC2, SC_III=SC3,
# SC_IV=SC1).
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/superclusters/SC_mapping.csv
#   data/superclusters/DEGs_TFs/SC<n>/SC<n>_gene_list.csv
#   data/ComPlEx/Orthogroups_20240613.tsv   (or data/annotation/…)
#   data/TranscriptionFactors/TF_<species>.RData   (for the TF-DEG pies)
# Outputs:
#   figures/Figure3C.png ; figures/Fig3C_pie_values.csv
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
NAMES <- c(A="Arabidopsis only", C="Coniferous only", D="Deciduous only",
           AC="Arabidopsis & Coniferous", AD="Arabidopsis & Deciduous",
           CD="Coniferous & Deciduous", ACD="All three groups")
pie_of <- function(orig, tfonly) {
  f <- Sys.glob(sprintf("data/superclusters/DEGs_TFs/%s/%s_gene_list.csv", orig, orig))
  sc <- sub("\\.\\d+$", "", read.csv(f[1], stringsAsFactors = FALSE)[, 2])
  if (tfonly) sc <- sc[sc %in% TF]
  ogi <- sapply(sc, function(g) { i <- mget(g, envir = g2og, ifnotfound = NA)[[1]]; if (is.null(i) || is.na(i)) NA else i })
  cl <- sapply(sc, clade); ok <- !is.na(ogi)
  ogset <- tapply(cl[ok], ogi[ok], function(v) paste(sort(unique(v)), collapse = ""))
  t <- table(factor(ogset, levels = CAT)); 100 * t / sum(t)
}

res <- list()
for (k in seq_len(nrow(mapping))) for (ty in c("DEG","TF-DEG")) {
  v <- pie_of(mapping$Original_supercluster[k], ty == "TF-DEG")
  res[[length(res)+1]] <- data.frame(SC = mapping$Main_text[k], type = ty,
                                     category = factor(NAMES[CAT], levels = NAMES[CAT]), pct = as.numeric(v))
}
D <- do.call(rbind, res)
dir.create("figures", showWarnings = FALSE)
write.csv(reshape(transform(D, category = names(NAMES)[match(category, NAMES)]),
                  idvar = c("SC","type"), timevar = "category", direction = "wide"),
          "figures/Fig3C_pie_values.csv", row.names = FALSE)

pal <- c("Arabidopsis only"="#3F936A","Coniferous only"="#786C9C","Deciduous only"="#E49C24",
         "Arabidopsis & Coniferous"="#6CB490","Arabidopsis & Deciduous"="#E48460",
         "Coniferous & Deciduous"="#9090B4","All three groups"="#A8C048")
D$SC <- factor(D$SC, levels = mapping$Main_text); D$type <- factor(D$type, levels = c("DEG","TF-DEG"))
p <- ggplot(D, aes(x = "", y = pct, fill = category)) +
  geom_col(width = 1, colour = "white", linewidth = 0.3) + coord_polar("y") +
  facet_grid(SC ~ type) + scale_fill_manual(values = pal, name = NULL) +
  theme_void(base_size = 12) + theme(legend.position = "right")
ggsave("figures/Figure3C.png", p, width = 7, height = 9, dpi = 120)
cat("Wrote figures/Figure3C.png and figures/Fig3C_pie_values.csv\n")
print(reshape(D, idvar = c("SC","type"), timevar = "category", direction = "wide"), row.names = FALSE)
