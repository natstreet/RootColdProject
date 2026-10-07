#!/usr/bin/env Rscript
# export_open_formats.R
#
# Re-export the R-serialised objects in this deposit as UTF-8, tab-separated
# text so that the resource is usable without R. Addresses FAIR I1
# (broadly applicable representation) and R1.3 (community formats).
#
# Run from the directory that contains data/ :
#     Rscript export_open_formats.R [outdir]
# Default outdir is data/open/
#
# Every output is written with row.names = FALSE and quote = FALSE. The
# deposit's own ComPlEx_MR_Network.R writes with col.names = TRUE and default
# row.names = TRUE, which emits a header one field short of its data rows and
# breaks naive parsers (pandas, csvkit, awk). Do not copy that pattern.

suppressPackageStartupMessages(library(tidyverse))
# arrow is needed only for step 5 (orthologs.parquet). Optional by design: the
# other seven exports must not fail because one reader is unavailable.
HAVE_ARROW <- requireNamespace("arrow", quietly = TRUE)

# The deposited ComPlEx_MR_Network.R refers to this file as data/orthologs.parquet
# (in a commented-out block); in the archive it sits at data/orthologs/orthologs.parquet.
# Accept either rather than assuming one. Resolved here because steps 4 and 5 both need it.
pq <- c("data/orthologs/orthologs.parquet", "data/orthologs.parquet")
pq <- pq[file.exists(pq)][1]

outdir <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(outdir)) outdir[1] else "data/open"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
wr <- function(x, name) {
  p <- file.path(outdir, name)
  write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("  %-42s %9d rows %6.1f MB\n", name, nrow(x), file.size(p) / 1e6))
  invisible(p)
}

SPECIES <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
# Species code used in the DE object filenames, mapped to the network/ComPlEx code
DE_CODE <- c(col0 = "col0", ost0 = "ost0", aspen = "pt", birch = "bp",
             spruce = "pa", pine = "ps")
TIMEPOINTS <- c("6h", "24h", "3d", "10d")

cat("1. conserved core (one gene per species per orthogroup)\n")
core <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")
clus <- ld("data/ComPlEx/cliques/clique_HMcluster_list.rda")
core_out <- core %>%
  left_join(as_tibble(clus), by = c("cliqueID" = "ID")) %>%
  relocate(OrthoGroup, cliqueID, Cluster)
wr(core_out, "conserved_core.tsv")

cat("2. co-expressolog pairs, all species pairs in one long table\n")
coex <- ld("data/ComPlEx/co_expressologs.RData")
wr(coex, "co_expressologs.tsv")

cat("3. per-pair comparison tables (neighbourhood overlap, effect size, p-values)\n")
# The Arabidopsis / Symbol / Name columns are a property of the orthogroup, not of the
# species pair, and repeating them across 1.35 M rows inflates this table roughly
# thirtyfold. They are dropped here and provided once in orthogroup_annotation.tsv;
# join on OrthoGroup to recover them.
pair_files <- list.files("data/ComPlEx/PairwiseComparisons", full.names = TRUE)
cmp <- map_dfr(pair_files, function(f) {
  t <- ld(f)
  sp <- gsub("^comparison-|-abspearsonMR0\\.03-table\\.RData$", "", basename(f))
  as_tibble(t) %>%
    dplyr::select(-any_of(c("Arabidopsis", "Symbol", "Name"))) %>%
    mutate(species_pair = sp, .before = 1)
})
wr(cmp, "pairwise_comparisons.tsv")

cat("4. orthogroup annotation (Arabidopsis locus, symbol, name) -- REBUILT, see note\n")
# DO NOT export data/ComPlEx/annotation-*.RData. Those tables, and the Arabidopsis /
# Symbol / Name columns of the PairwiseComparisons tables, were built from
# data/annotation/Orthogroups_130323_predefined_tree.tsv, whereas the co-expression
# analysis itself used the orthology in orthologs.parquet. Both orthology runs number
# their groups OG0000000, OG0000001, ... so the identifiers collide while meaning
# different gene sets, and the annotation is therefore attached to the wrong groups.
# Measured: the deposited annotation's Arabidopsis locus is a member of the
# same-numbered analysis orthogroup in 2 of 6566 cases (1 of 95 conserved-core groups).
# Worked example: OG0000065 is labelled CBF1/CBF2/CBF3/DREB1 by the deposited
# annotation but actually contains the WRKY13/23/24/28/43/45 group -- which would lead
# a reader to contradict the paper's finding that ICE1-CBF-COR is absent from the
# conserved core.
#
# This step rebuilds the mapping from the analysis orthology plus the TAIR alias table,
# and asserts that every annotated locus really is a member of the group it annotates.
if (!HAVE_ARROW || is.na(pq)) {
  cat("  SKIPPED: needs arrow + orthologs.parquet to rebuild correctly.\n")
  cat("  Do NOT substitute data/ComPlEx/annotation-*.RData: it is mis-joined.\n")
} else {
  ath <- arrow::read_parquet(pq) %>%
    filter(species == "Arabidopsis_thaliana") %>%
    transmute(OrthoGroup = Ortholog_Group,
              locus = toupper(trimws(as.character(gene)))) %>%
    filter(grepl("^AT[1-5CM]G[0-9]{5}$", locus)) %>%
    distinct()
  alias <- read.delim("data/annotation/gene_aliases_20140331.txt",
                      stringsAsFactors = FALSE, colClasses = "character") %>%
    transmute(locus = toupper(trimws(locus_name)), symbol, full_name) %>%
    group_by(locus) %>%
    summarise(Symbol = paste(sort(unique(symbol[nzchar(symbol)])), collapse = "; "),
              Name   = paste(sort(unique(full_name[nzchar(full_name)])), collapse = "; "),
              .groups = "drop")
  annot <- ath %>%
    left_join(alias, by = "locus") %>%
    group_by(OrthoGroup) %>%
    summarise(Arabidopsis = paste(sort(unique(locus)), collapse = "; "),
              Symbol = substr(paste(sort(unique(Symbol[!is.na(Symbol) & nzchar(Symbol)])),
                                    collapse = "; "), 1, 300),
              Name   = substr(paste(sort(unique(Name[!is.na(Name) & nzchar(Name)])),
                                    collapse = "; "), 1, 300),
              .groups = "drop")
  # self-consistency: the first annotated locus must belong to the group it annotates
  memb <- split(ath$locus, ath$OrthoGroup)
  first_at <- sub(";.*$", "", annot$Arabidopsis)
  ok <- mapply(function(og, l) l %in% memb[[og]], annot$OrthoGroup, first_at)
  stopifnot(all(ok))
  cat(sprintf("  self-consistency check: %d/%d annotated loci are members\n", sum(ok), length(ok)))
  wr(annot, "orthogroup_annotation.tsv")
}

cat("5. gene <-> orthogroup map, all species\n")
if (!HAVE_ARROW) {
  cat("  SKIPPED: install.packages('arrow') to export orthologs.parquet\n")
} else if (is.na(pq)) {
  cat("  SKIPPED: orthologs.parquet not found at either documented path\n")
} else {
  wr(as_tibble(arrow::read_parquet(pq)), "gene_orthogroup_map.tsv")
}

cat("6. conserved vs species-specific orthogroup gene sets\n")
sets <- map_dfr(SPECIES, function(s) {
  f_common <- sprintf("data/ComPlEx/Common_OG_genes_%s_005.RData", s)
  f_unique <- sprintf("data/ComPlEx/Unique_OG_genes_%s_005.RData", s)
  bind_rows(
    if (file.exists(f_common)) tibble(species = s, class = "conserved", gene = ld(f_common)),
    if (file.exists(f_unique)) tibble(species = s, class = "species_specific", gene = ld(f_unique))
  )
})
wr(sets, "conserved_vs_specific_genes.tsv")

cat("7. differential expression, long format\n")
de <- map_dfr(SPECIES, function(s) {
  cd <- DE_CODE[[s]]
  f_l <- sprintf("data/DEGs/log2FoldChange_%s.RData", cd)
  f_p <- sprintf("data/DEGs/padj_%s.RData", cd)
  if (!file.exists(f_l) || !file.exists(f_p)) return(NULL)
  L <- as.data.frame(ld(f_l)); P <- as.data.frame(ld(f_p))
  g <- intersect(rownames(L), rownames(P))
  map_dfr(seq_along(TIMEPOINTS), function(i) {
    lc <- grep(TIMEPOINTS[i], colnames(L), fixed = TRUE)[1]
    pc <- grep(TIMEPOINTS[i], colnames(P), fixed = TRUE)[1]
    if (is.na(lc) || is.na(pc)) return(NULL)
    tibble(species = s, gene = g, timepoint = TIMEPOINTS[i],
           log2FoldChange = L[g, lc], padj = P[g, pc])
  })
})
wr(de, "differential_expression.tsv")

cat("8. transcription factors (schemas harmonised, see note)\n")
# The five deposited TF objects do not share a schema: TF_col0.RData has named columns
# (TF_ID, Gene_ID, Family) while TF_aspen/birch/spruce/pine.RData have unnamed columns
# V1 (gene) and V2 (family). Row-binding them without harmonising silently puts the
# conifer and Populus gene identifiers in a different column from the Arabidopsis ones,
# which makes any join against them return zero rows. Harmonise to species/gene/family.
# Note also that there is no TF_ost0 object: Ost-0 is Arabidopsis thaliana and shares
# the Col-0 transcription-factor set.
tf <- map_dfr(c("col0", "aspen", "birch", "spruce", "pine"), function(s) {
  f <- sprintf("data/TranscriptionFactors/TF_%s.RData", s)
  if (!file.exists(f)) return(NULL)
  o <- as.data.frame(ld(f), stringsAsFactors = FALSE)
  if (all(c("Gene_ID", "Family") %in% names(o))) {
    tibble(species = s, gene = as.character(o$Gene_ID), family = as.character(o$Family))
  } else {
    tibble(species = s, gene = as.character(o[[1]]), family = as.character(o[[2]]))
  }
}) %>% filter(!is.na(gene), nzchar(gene)) %>% distinct()
# Ost-0 shares the Col-0 TF set; make that explicit rather than leaving a gap.
tf <- bind_rows(tf, tf %>% filter(species == "col0") %>% mutate(species = "ost0"))
wr(tf, "transcription_factors.tsv")

cat("\nwrote to ", normalizePath(outdir), "\n", sep = "")
