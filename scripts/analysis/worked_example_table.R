#!/usr/bin/env Rscript
# worked_example_table.R — the Fig. S1 worked example in numbers: per genotype, the clique member, its Methods DEG call
# and direction (tables/core_response_patterns.tsv), its number of conserved neighbours, how many are positively and
# negatively correlated with it (tables/core_conserved_neighbours.tsv.gz), and the predicted transcription factors
# among them (data/TranscriptionFactors/TF_<species>.RData; Ost-0 uses the Arabidopsis list), with their families.
# Usage (repository root): Rscript scripts/analysis/worked_example_table.R <orthogroup> <out tsv>
suppressMessages(library(data.table))
gzread <- function(f) fread(cmd = paste("gzip -dc", shQuote(f)))   # gzipped tables without R.utils
a <- commandArgs(trailingOnly = TRUE); OG <- a[1]
SP <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
strip <- function(x) sub("[.][0-9]+$", "", as.character(x))
pat <- fread("tables/core_response_patterns.tsv")[OrthoGroup == OG]
cn <- gzread("tables/core_conserved_neighbours.tsv.gz")[OrthoGroup == OG]
if (nrow(pat) != 1 || !nrow(cn)) stop(OG, " is not in the core tables")
tf <- list()
for (s in c("col0", "aspen", "birch", "spruce", "pine")) { x <- as.data.frame(ld(sprintf("data/TranscriptionFactors/TF_%s.RData", s)))
  fam <- if (all(c("Gene_ID", "Family") %in% names(x))) setNames(as.character(x$Family), strip(x$Gene_ID)) else setNames(as.character(x[[2]]), strip(x[[1]]))
  tf[[s]] <- fam[!duplicated(names(fam))] }
tf$ost0 <- tf$col0
out <- rbindlist(lapply(SP, function(s) {
  c_s <- cn[species == s]; f <- tf[[s]][strip(c_s$neighbour)]; isTF <- !is.na(f)
  data.table(genotype = s, member = c_s$member[1], response = pat[[paste0("response_", s)]],
             timepoints_DE = pat[[paste0("timepoints_DE_", s)]], conserved_neighbours = nrow(c_s),
             positive = sum(c_s$r > 0), negative = sum(c_s$r < 0), transcription_factors = sum(isTF),
             tf_neighbours = paste(sprintf("%s (%s)", c_s$neighbour[isTF], f[isTF]), collapse = "; "))
}))
fwrite(out, a[2], sep = "\t"); print(out[, -"tf_neighbours"])
