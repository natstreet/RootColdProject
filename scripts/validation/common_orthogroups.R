#!/usr/bin/env Rscript
# common_orthogroups.R — the commonly differentially expressed orthogroups (a DEG in all six genotypes): how many
# there are under each gene-to-orthogroup map, and how large they are.
#
# The common set is taken from data/DEGs/og_summary.RData (n_species_with_DEGs == 6), which is what Table S1 uses.
# Its reference set is the orthogroups with genes in all six genotypes in the same table (n_species_with_genes == 6).
# Orthogroup size is the number of genes over the five species in Orthogroups_20240613.tsv (Col-0 and Ost-0 share
# the Arabidopsis column). The same common-set count is also computed by mapping each genotype's DEGs through the
# other maps in the repository and the deposit, to show which map the published set corresponds to.
#
# Inputs (unpack data.tar.gz from figshare into ./data):
#   data/DEGs/og_summary.RData, data/DEGs/DEGs_<genotype>.RData
#   data/*/Orthogroups_20240613.tsv
#   data/ComPlEx/gene2OG/gene2OG_<genotype>.rda                 (in this repository)
#   data/open/gene_orthogroup_map.tsv                           (deposit v4 or later; skipped if absent)
# Output: results/common_orthogroups.tsv (quantity, value)
#
# Usage: Rscript scripts/validation/common_orthogroups.R

ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
GT <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
SPCOL <- c(col0 = "Arabidopsis_thaliana", ost0 = "Arabidopsis_thaliana", aspen = "Populus_tremula",
           birch = "Betula_pendula", spruce = "Picea_abies", pine = "Pinus_sylvestris")
MIN_SIZE <- 20

og <- read.delim(Sys.glob("data/*/Orthogroups_20240613.tsv")[1], check.names = FALSE, stringsAsFactors = FALSE)
members <- lapply(unique(SPCOL), function(s) { x <- og[[s]]; x[is.na(x)] <- ""; strsplit(x, ",\\s*") })
names(members) <- unique(SPCOL)
size <- setNames(Reduce(`+`, lapply(members, lengths)), og$Ortholog_Group)

os <- ld("data/DEGs/og_summary.RData")
common <- os$OrthoGroup[os$n_species_with_DEGs == 6]
ref <- os$OrthoGroup[os$n_species_with_genes == 6]
stopifnot(all(common %in% ref), all(ref %in% names(size)))

out <- data.frame(
  quantity = c("common_orthogroups", "common_with_at_least_20_genes", "common_pct_at_least_20_genes",
               "orthogroups_with_genes_in_all_six_genotypes", "of_those_with_at_least_20_genes",
               "of_those_pct_at_least_20_genes"),
  value = c(length(common), sum(size[common] >= MIN_SIZE), 100 * mean(size[common] >= MIN_SIZE),
            length(ref), sum(size[ref] >= MIN_SIZE), 100 * mean(size[ref] >= MIN_SIZE)))

# the common-set count under each gene -> orthogroup map
deg <- lapply(setNames(GT, GT), function(g) unique(as.character(ld(sprintf("data/DEGs/DEGs_%s.RData", g)))))
count_common <- function(g2og) {                   # g2og: genotype -> named vector gene -> orthogroup
  sets <- lapply(GT, function(g) unique(na.omit(g2og[[g]][deg[[g]]])))
  length(Reduce(intersect, sets))
}
tsv_map <- lapply(setNames(GT, GT), function(g) {
  m <- members[[SPCOL[[g]]]]; setNames(rep(og$Ortholog_Group, lengths(m)), unlist(m)) })
out <- rbind(out, data.frame(quantity = "common_via_Orthogroups_20240613.tsv", value = count_common(tsv_map)))

g2f <- Sys.glob("data/ComPlEx/gene2OG/gene2OG_*.rda")
if (length(g2f) == 6) {
  g2 <- lapply(setNames(GT, GT), function(g) {
    x <- as.data.frame(ld(g2f[tolower(sub("^gene2OG_(.*)[.]rda$", "\\1", basename(g2f))) == g]))
    setNames(as.character(x$OrthoGroup), as.character(x[[2]])) })       # columns: OrthoGroup, <genotype> gene
  out <- rbind(out, data.frame(quantity = "common_via_ComPlEx_gene2OG", value = count_common(g2)))
}
gm <- "data/open/gene_orthogroup_map.tsv"
if (file.exists(gm)) {
  m <- read.delim(gm, stringsAsFactors = FALSE)
  gmap <- lapply(setNames(GT, GT), function(g) { x <- m[m$species == SPCOL[[g]], ]; setNames(x$Ortholog_Group, x$gene) })
  out <- rbind(out, data.frame(quantity = "common_via_gene_orthogroup_map.tsv", value = count_common(gmap)))
}

dir.create("results", showWarnings = FALSE)
write.table(out, "results/common_orthogroups.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
