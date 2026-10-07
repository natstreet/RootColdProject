# tf_cliques.R — the transcription-factor cliques of the conserved core (Figure 5, tf_clique_directions.R): the
# core cliques in which at least MIN_TF of the six members are predicted transcription factors
# (data/TranscriptionFactors/TF_<species>.RData; Ost-0 uses the Arabidopsis list). Five of six is the rule that
# reproduces the published selection (one bZIP61 clique member, in birch, is not in the birch prediction).
# Each clique is labelled by its Col-0 member's symbol (gene_aliases_20140331.txt; a leading "AT" or "AT-" is
# dropped and BZIP is written bZIP), with the Ost-0 member's symbol added ("ANAC018/ANAC002") where the two differ;
# a gene without a symbol is labelled by its locus. Cliques are returned in alphabetical order of label.
# Usage: source("scripts/analysis/tf_cliques.R"); tfc <- tf_cliques()   (columns tf, OrthoGroup, col0 ... pine, n_tf)
tf_cliques <- function(MIN_TF = 5) {
  sp <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
  strip <- function(x) sub("[.][0-9]+$", "", as.character(x))
  core <- as.data.frame(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
  tf <- list()
  for (s in c("col0", "aspen", "birch", "spruce", "pine")) {
    e <- new.env(); load(sprintf("data/TranscriptionFactors/TF_%s.RData", s), envir = e); x <- get(ls(e)[1], e)
    tf[[s]] <- unique(c(strip(x[[1]]), strip(x[[2]])))
  }
  tf$ost0 <- tf$col0
  n_tf <- rowSums(sapply(sp, function(s) strip(core[[s]]) %in% tf[[s]]))
  al <- read.delim("data/annotation/gene_aliases_20140331.txt", colClasses = "character")
  sym <- function(g) {
    v <- al$symbol[toupper(al$locus_name) == strip(g)]
    if (!length(v) || !nzchar(v[1])) return(strip(g))
    s <- sub("^AT-?(?=[A-Za-z])", "", v[1], perl = TRUE)
    sub("^BZIP", "bZIP", s)
  }
  r <- core[n_tf >= MIN_TF, c("OrthoGroup", sp)]
  if (!nrow(r)) stop("no core clique has at least ", MIN_TF, " transcription-factor members")
  r$tf <- mapply(function(c0, o0) if (strip(c0) == strip(o0)) sym(c0) else paste0(sym(c0), "/", sym(o0)), r$col0, r$ost0)
  r$n_tf <- n_tf[n_tf >= MIN_TF]
  r <- r[order(toupper(r$tf)), c("tf", "OrthoGroup", sp, "n_tf")]
  rownames(r) <- NULL
  r
}
