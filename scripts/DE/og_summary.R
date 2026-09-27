#!/usr/bin/env Rscript
# og_summary.R — per orthogroup, the genotypes with a gene in it and the genotypes with a DEG in it.
#
# Inputs:
#   data/ComPlEx/gene2OG/gene2OG_<Genotype>.rda   gene -> orthogroup, per genotype (this repository)
#   data/DEGs/DEGs_<genotype>.RData               the DEGs of each genotype (union over timepoints; deposit)
# Output: data/DEGs/og_summary.RData (object og_summary: OrthoGroup, species_with_genes, species_with_DEGs (lists, in
#         the genotype order col0, ost0, aspen, birch, spruce, pine), n_species_with_genes, n_species_with_DEGs), or
#         the path given as the first argument.
# Read by: make_Fig1.R, make_TableS1.R, common_orthogroups.R
# Usage:  Rscript scripts/DE/og_summary.R [output.RData]
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else "data/DEGs/og_summary.RData"
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

sp <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
fn <- c(col0 = "Col0", ost0 = "Ost0", aspen = "aspen", birch = "birch", spruce = "spruce", pine = "pine")
M <- bind_rows(lapply(sp, function(s) {
  d <- ld(sprintf("data/ComPlEx/gene2OG/gene2OG_%s.rda", fn[s])); names(d)[1:2] <- c("OrthoGroup", "gene")
  data.frame(OrthoGroup = d$OrthoGroup, gene = d$gene, species = s) }))
degs <- setNames(lapply(sp, function(s) {
  v <- ld(sprintf("data/DEGs/DEGs_%s.RData", s))
  if (is.list(v) && !is.data.frame(v)) v <- unlist(v)
  if (is.data.frame(v)) v <- v[[1]]
  unique(as.character(v)) }), sp)
M$deg <- mapply(function(g, s) g %in% degs[[s]], M$gene, M$species)
og_summary <- M %>% group_by(OrthoGroup) %>%
  summarise(species_with_genes = list(sp[sp %in% unique(species)]),
            species_with_DEGs  = list(sp[sp %in% unique(species[deg])]), .groups = "drop") %>%
  mutate(n_species_with_genes = lengths(species_with_genes), n_species_with_DEGs = lengths(species_with_DEGs))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
save(og_summary, file = out)
cat("og_summary:", nrow(og_summary), "orthogroups ->", out, "\n")
