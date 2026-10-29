#!/usr/bin/env Rscript
# core_status_columns.R — deposit, version 4 (ledger R037 item 2, R042 item 1; DECISIONS.md D050, D066): each orthogroup's
# place relative to the conserved core, the sample-stability label of the core orthogroups, and the near-core list.
#   core_status  core (counted as conserved: complete clique and the D026 size rule) | set aside by size (complete clique,
#                not counted) | near-core (co-expressolog in all 15 comparisons, Max.p.val < 0.1, but no complete
#                clique) | none
#   stability    stable | sample-sensitive (core orthogroups only; D066, from core_stability_label.R) ; share_passing =
#                passing / complete cliques under the leave-one-sample-out test
# Inputs (working tree): results/orthogroup_size_rule.tsv (eligible orthogroups, expressed_genes; D080), results/orthogroup_sizes_all.tsv (all15), results/core_stability_label.tsv,
#   data/ComPlEx/PairwiseComparisons/comparison-col0-*-table.RData (Col-0 genes tested per orthogroup and their symbols)
# Outputs (<open dir>): orthogroup_core_status.tsv (every eligible orthogroup), near_core_orthogroups.tsv; and the columns
#   core_status, stability, share_passing added to orthogroup_annotation.tsv and conserved_core*.tsv when present
# Usage (working tree): Rscript scripts/deposit/core_status_columns.R <open dir>
suppressPackageStartupMessages(library(data.table))
O <- commandArgs(trailingOnly = TRUE)[1]
ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(get(ls(e)[1], e)) }
# fwrite writes an empty string as "" (to tell it from NA); the deposit's tables leave empty fields empty, so empty strings
# become NA before writing with na = "" (ledger R046)
blank_na <- function(t) { for (k in names(t)) if (is.character(t[[k]])) set(t, which(t[[k]] == ""), k, NA_character_); t }
sz <- fread("results/orthogroup_sizes_all.tsv"); sr <- fread("results/orthogroup_size_rule.tsv"); lab <- fread("results/core_stability_label.tsv")
# eligible orthogroups and expressed_genes as the size rule defines them (D080: genes passing the ComPlEx expression filter, at least
# one per species); all15 from the size-flag table, whose eligible set contains every D080-eligible orthogroup (checked below)
x <- merge(sr[, .(OrthoGroup, expressed_genes, has_complete_clique, counted_as_conserved)], sz[, .(OrthoGroup, all15)], by = "OrthoGroup", all.x = TRUE)
if (anyNA(x$all15)) stop("D080-eligible orthogroups missing from orthogroup_sizes_all.tsv")
x[, core_status := fifelse(counted_as_conserved, "core", fifelse(has_complete_clique, "set aside by size", fifelse(all15, "near-core", "none")))]
x <- merge(x, lab[, .(OrthoGroup, stability = core_status, share_passing)], by = "OrthoGroup", all.x = TRUE)
if (x[core_status == "core", .N] != nrow(lab) || x[core_status == "core", anyNA(stability)]) stop("core orthogroups and stability labels do not agree")
fwrite(x[order(OrthoGroup)], file.path(O, "orthogroup_core_status.tsv"), sep = "\t", na = "")
print(x[, .N, by = .(core_status, stability)][order(core_status)])
nc <- x[core_status == "near-core", OrthoGroup]
ann <- unique(rbindlist(lapply(list.files("data/ComPlEx/PairwiseComparisons", "^comparison-col0-.*-table\\.RData$", full.names = TRUE),
  function(f) ld(f)[OrthoGroup %in% nc, .(OrthoGroup, gene = Species1, symbol = Symbol)])))
ann[is.na(symbol) | symbol == "", symbol := NA_character_]
nl <- ann[, .(col0_genes_tested = uniqueN(gene), col0_symbols = paste(sort(unique(na.omit(symbol))), collapse = " | ")), by = OrthoGroup]
nl <- merge(x[core_status == "near-core", .(OrthoGroup, expressed_genes)], nl, by = "OrthoGroup", all.x = TRUE)
fwrite(blank_na(nl[order(OrthoGroup)]), file.path(O, "near_core_orthogroups.tsv"), sep = "\t", na = "")
cat("near-core orthogroups:", nrow(nl), "\n")
for (f in c("orthogroup_annotation.tsv", list.files(O, "^conserved_core.*\\.tsv$"))) { p <- file.path(O, f); if (!file.exists(p)) next
  # every field as text, unstripped and unquoted, so the existing columns are written back unchanged (ledger R046)
  t <- fread(p, sep = "\t", colClasses = "character", strip.white = FALSE, quote = "", na.strings = NULL)
  k <- intersect(c("OrthoGroup", "orthogroup", "Orthogroup"), names(t))[1]; if (is.na(k)) { cat("no orthogroup column in", f, "\n"); next }
  t[, intersect(c("core_status", "stability", "share_passing"), names(t)) := NULL]
  t <- merge(t, x[, .(OrthoGroup, core_status, stability, share_passing)], by.x = k, by.y = "OrthoGroup", all.x = TRUE, sort = FALSE)
  t[is.na(core_status), core_status := "none"]; fwrite(blank_na(t), p, sep = "\t", na = ""); cat(f, ": columns added\n") }
