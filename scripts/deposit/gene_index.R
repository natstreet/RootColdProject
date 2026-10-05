#!/usr/bin/env Rscript
# gene_index.R — the deposit's gene_index.tsv: one row per gene per genotype, joining the other open tables so a
# single gene can be looked up without a join (DATA_DICTIONARY.md). Until now the table had no producer; this one
# rebuilds it from the open tables in the same directory.
#
# Rows: every (species, gene) in differential_expression.tsv (one row per gene; four timepoints each).
# Columns, and the table each comes from:
#   orthogroup                            gene_orthogroup_map.tsv (species name -> genotype; Arabidopsis for col0 and ost0)
#   in_conserved_core, core_cluster, core_cluster_manuscript, core_cliqueID, core_conservedType
#                                         conserved_core.tsv (the gene is the genotype's member of a core clique)
#   conservation_class                    conserved_vs_specific_genes.tsv
#   n_co_expressologs                     co_expressologs.tsv (pairs in which the gene takes part, either side)
#   is_TF, tf_family                      transcription_factors.tsv (families sorted and joined with "; ")
#   arabidopsis_locus/_symbol/_name       orthogroup_annotation.tsv, by orthogroup
#   responds_to_cold, n_timepoints_significant, peak_timepoint, peak_lfc, max_abs_lfc, lfc_*, padj_*, padj_lfc0_*
#                                         differential_expression.tsv (is_deg is the article's call; the peak is the
#                                         timepoint of largest |log2FoldChange|, whatever its significance;
#                                         peak_lfc and max_abs_lfc rounded to 3 decimals, as deposited)
#
# Usage: Rscript scripts/deposit/gene_index.R <data/open dir> [output file, default <dir>/gene_index.tsv]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
open <- a[1]; stopifnot(!is.na(open), dir.exists(open))
out <- if (length(a) > 1) a[2] else file.path(open, "gene_index.tsv")
rd <- function(f) fread(file.path(open, f), sep = "\t", colClasses = "character", na.strings = NULL, quote = "")
TP <- c("6h", "24h", "3d", "10d")
SPECIES <- c(col0 = "Arabidopsis_thaliana", ost0 = "Arabidopsis_thaliana", aspen = "Populus_tremula",
             birch = "Betula_pendula", spruce = "Picea_abies", pine = "Pinus_sylvestris")
num <- function(x) as.numeric(ifelse(x %in% c("", "NA"), NA, x))

de <- rd("differential_expression.tsv")
de[, lfc := num(log2FoldChange)]
if (!all(de$timepoint %in% TP)) stop("unexpected timepoints")
g <- unique(de[, .(species, gene)])

map <- rd("gene_orthogroup_map.tsv")[species %in% SPECIES, .(sp_name = species, gene, orthogroup = Ortholog_Group)]
map <- unique(map)
if (anyDuplicated(map[, .(sp_name, gene)])) stop("a gene has more than one orthogroup")
g[, sp_name := SPECIES[species]]
g <- merge(g, map, by = c("sp_name", "gene"), all.x = TRUE, sort = FALSE)

core <- rd("conserved_core.tsv")
cm <- rbindlist(lapply(names(SPECIES), function(s)
  core[, .(species = s, gene = get(s), core_cluster = Cluster,
           core_cluster_manuscript = if ("ManuscriptCluster" %in% names(core)) ManuscriptCluster else "",
           core_cliqueID = cliqueID, core_conservedType = conservedType)]))
if (anyDuplicated(cm[, .(species, gene)])) stop("a gene is the member of two core cliques")
g <- merge(g, cm, by = c("species", "gene"), all.x = TRUE, sort = FALSE)
g[, in_conserved_core := !is.na(core_cliqueID)]

cs <- rd("conserved_vs_specific_genes.tsv")
if (anyDuplicated(cs[, .(species, gene)])) stop("a gene has two conservation classes")
g <- merge(g, cs[, .(species, gene, conservation_class = class)], by = c("species", "gene"), all.x = TRUE, sort = FALSE)

cx <- rd("co_expressologs.tsv")
nco <- rbind(cx[, .(species = Species1, gene = GeneSpecies1)], cx[, .(species = Species2, gene = GeneSpecies2)])[, .N, by = .(species, gene)]
g <- merge(g, nco[, .(species, gene, n_co_expressologs = N)], by = c("species", "gene"), all.x = TRUE, sort = FALSE)
g[is.na(n_co_expressologs), n_co_expressologs := 0L]

tf <- rd("transcription_factors.tsv")[, .(tf_family = paste(sort(unique(family)), collapse = "; ")), by = .(species, gene)]
g <- merge(g, tf, by = c("species", "gene"), all.x = TRUE, sort = FALSE)
g[, is_TF := !is.na(tf_family)]

an <- rd("orthogroup_annotation.tsv")[, .(orthogroup = OrthoGroup, arabidopsis_locus = Arabidopsis,
                                          arabidopsis_symbol = Symbol, arabidopsis_name = Name)]
g <- merge(g, an, by = "orthogroup", all.x = TRUE, sort = FALSE)

de[, deg := toupper(is_deg) == "TRUE"]
s1 <- de[, .(responds_to_cold = any(deg), n_timepoints_significant = sum(deg),
             peak_timepoint = if (all(is.na(lfc))) NA_character_ else timepoint[which.max(abs(lfc))],
             peak_lfc = if (all(is.na(lfc))) NA_real_ else lfc[which.max(abs(lfc))]), by = .(species, gene)]
s1[, peak_lfc := round(peak_lfc, 3)][, max_abs_lfc := abs(peak_lfc)]   # 3 decimals, as deposited
w <- dcast(de, species + gene ~ timepoint, value.var = c("log2FoldChange", "padj", "padj_lfc0"))
setnames(w, paste0("log2FoldChange_", TP), paste0("lfc_", TP))
setnames(w, paste0("padj_lfc0_", TP), paste0("padj_lfc0_", TP))
g <- merge(g, s1, by = c("species", "gene"), all.x = TRUE, sort = FALSE)
g <- merge(g, w, by = c("species", "gene"), all.x = TRUE, sort = FALSE)

cols <- c("species", "gene", "orthogroup", "in_conserved_core", "core_cluster", "core_cluster_manuscript",
          "core_cliqueID", "core_conservedType", "conservation_class", "n_co_expressologs", "is_TF", "tf_family",
          "arabidopsis_locus", "arabidopsis_symbol", "arabidopsis_name", "responds_to_cold",
          "n_timepoints_significant", "peak_timepoint", "peak_lfc", "max_abs_lfc",
          paste0("lfc_", TP), paste0("padj_", TP), paste0("padj_lfc0_", TP))
g <- g[, ..cols]
setorder(g, species, gene)                                     # rows sorted by genotype, then gene
for (b in c("in_conserved_core", "is_TF", "responds_to_cold")) set(g, j = b, value = ifelse(g[[b]], "True", "False"))
fwrite(g, out, sep = "\t", quote = FALSE, na = "")
cat(sprintf("gene_index.tsv: %d rows, %d columns; core members per genotype: %s -> %s\n", nrow(g), ncol(g),
            paste(g[in_conserved_core == "True", .N, by = species][, paste0(species, " ", N)], collapse = ", "), out))
