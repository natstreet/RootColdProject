#!/usr/bin/env Rscript
# cluster_profiles.R — cluster-level cold response of the conserved-core cliques (Results, clusters 1 and 2).
#
# Produces the numbers quoted in the Results for manuscript clusters 1 and 2: per cluster, genotype and timepoint, the
# median log2 fold change of the clique member genes, the equivalent fold, and how many genes rise, fall and change at
# least two-fold. Also lists cluster 2's members with their Arabidopsis annotation, which is the basis for naming the
# member classes in the text.
#
# Inputs (deposit, data/open/): gene_index.tsv, conserved_core.tsv, orthogroup_annotation.tsv.
# Cluster numbers: conserved_core.tsv's Cluster is the analysis numbering; the manuscript numbering is MAP in
# scripts/GO/make_TableS4.R, read from there.
# Fold changes are the deposit's log2FoldChange (DESeq2 MLE, not shrunken).
# Outputs: <out>/cluster_profiles.tsv, <out>/cluster2_members.tsv
#
# Usage: Rscript scripts/analysis/cluster_profiles.R <deposit data/open dir> <output dir>

a <- commandArgs(trailingOnly = TRUE); OPEN <- a[1]; OUT <- a[2]
stopifnot(!is.na(OPEN), !is.na(OUT)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SP <- c("col0", "ost0", "aspen", "birch", "spruce", "pine"); TPS <- c("6h", "24h", "3d", "10d")
src <- readLines("scripts/GO/make_TableS4.R")
MAP <- eval(parse(text = sub("^MAP <- ", "", regmatches(src, regexpr("MAP <- c\\(.*\\)", src)))))  # manuscript -> analysis
rd <- function(f) read.delim(file.path(OPEN, f), stringsAsFactors = FALSE, quote = "", check.names = FALSE)

gi <- rd("gene_index.tsv"); core <- rd("conserved_core.tsv"); ann <- rd("orthogroup_annotation.tsv")
core$ms_cluster <- as.integer(names(MAP))[match(core$Cluster, MAP)]
stopifnot(!anyNA(core$ms_cluster))
cat("orthogroups per manuscript cluster:", paste(sprintf("%s=%d", names(table(core$ms_cluster)), table(core$ms_cluster)),
                                                collapse = " "), "\n")

# clique member genes: exactly the gene named per genotype in conserved_core.tsv
long <- do.call(rbind, lapply(SP, function(s) data.frame(OrthoGroup = core$OrthoGroup, ms_cluster = core$ms_cluster,
                                                          species = s, gene = core[[s]], stringsAsFactors = FALSE)))
long <- long[!is.na(long$gene) & nzchar(long$gene), ]
key <- function(s, g) paste(s, g)
stopifnot(!anyDuplicated(key(gi$species, gi$gene)))
m <- cbind(long, gi[match(key(long$species, long$gene), key(gi$species, gi$gene)), paste0("lfc_", TPS)])
stopifnot(mean(rowSums(is.na(m[paste0("lfc_", TPS)])) == 0) > 0.99)

rows <- list()
for (c in sort(unique(m$ms_cluster))) for (sp in SP) {
  d <- m[m$ms_cluster == c & m$species == sp, ]
  if (!nrow(d)) next
  for (t in TPS) {
    x <- d[[paste0("lfc_", t)]]; x <- x[!is.na(x)]; md <- median(x)
    rows[[length(rows) + 1]] <- data.frame(
      ms_cluster = c, species = sp, timepoint = t, n = length(x), median_log2fc = round(md, 3),
      fold = if (md >= 0) round(2^md, 2) else round(-2^-md, 2),
      n_up = sum(x > 0), n_down = sum(x < 0), n_abs_log2fc_ge1 = sum(abs(x) >= 1))
  }
}
P <- do.call(rbind, rows)
write.table(P, file.path(OUT, "cluster_profiles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mem <- core[core$ms_cluster == 2, c("OrthoGroup", SP)]
mem <- merge(mem, ann[, c("OrthoGroup", "Arabidopsis", "Symbol", "Name")], by = "OrthoGroup", all.x = TRUE, sort = FALSE)
mem <- mem[match(core$OrthoGroup[core$ms_cluster == 2], mem$OrthoGroup), ]
write.table(mem, file.path(OUT, "cluster2_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("cluster_profiles.tsv: %d rows | cluster2_members.tsv: %d orthogroups\n", nrow(P), nrow(mem)))
