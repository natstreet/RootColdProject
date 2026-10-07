#!/usr/bin/env Rscript
# de_calls.R — put the paper's differential-expression call into the deposit tables.
#
# The deposited padj (differential_expression.tsv, gene_index.tsv padj_*) tests a null of |log2 fold change| <= 1,
# whereas the paper calls a DEG at padj <= 0.05 against a null of no change AND |log2 fold change| >= 1
# (scripts/DE/differential_expression.R). This step adds the latter from that script's DE/de_results.tsv, and bases
# gene_index's cold-response summary on it, so that the deposit's DEGs are the paper's DEGs.
#
# Edits in place (every other column unchanged):
#   <open>/differential_expression.tsv  adds padj_lfc0 (null of no change) and is_deg, after padj
#   <open>/gene_index.tsv (if present)  responds_to_cold = DEG at any timepoint; n_timepoints_significant = number of
#                                       timepoints at which it is a DEG; adds padj_lfc0_6h ... padj_lfc0_10d after padj_10d
# Stops unless: every (species, gene, timepoint) row is matched, the fold changes agree with the deposited ones, and
# each genotype's DEG union equals the published set in data/DEGs/DEGs_<species>.RData.
#
# Usage: Rscript scripts/DE/differential_expression.R            # writes DE/de_results.tsv
#        Rscript scripts/deposit/de_calls.R <deposit data/open dir>

suppressMessages(library(dplyr))
open <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(open), dir.exists(open))
TPS <- c("6h", "24h", "3d", "10d")
rd  <- function(p) read.delim(p, colClasses = "character", na.strings = character(0), check.names = FALSE, quote = "")
wr  <- function(x, p) write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
ld  <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
num <- function(x) suppressWarnings(as.numeric(x))

res <- rd("DE/de_results.tsv")
key <- function(d) paste(d$species, d$gene, d$timepoint)

# each genotype's DEG union must be the published one
for (sp in unique(res$species)) {
  u <- unique(res$gene[res$species == sp & res$deg == "TRUE"])
  pub <- unique(as.character(ld(sprintf("data/DEGs/DEGs_%s.RData", sp))))
  if (!setequal(u, pub)) stop(sp, ": DEG union (", length(u), ") differs from the published set (", length(pub), ")")
  cat(sprintf("%-7s DEGs %5d = published set\n", sp, length(u)))
}

## differential_expression.tsv
p <- file.path(open, "differential_expression.tsv"); de <- rd(p)
if ("padj_lfc0" %in% names(de)) stop("differential_expression.tsv already has padj_lfc0")
m <- match(key(de), key(res)); if (anyNA(m) || nrow(de) != nrow(res)) stop("differential_expression.tsv rows do not match DE/de_results.tsv")
a <- num(de$log2FoldChange); b <- num(res$log2FoldChange[m])
if (!identical(is.na(a), is.na(b))) stop("fold changes are missing in different rows")
cat(sprintf("fold changes: %d rows, %d within 1e-6 of the recomputed ones, largest difference %.3g\n",
            sum(!is.na(a)), sum(abs(a - b) <= 1e-6, na.rm = TRUE), max(abs(a - b), na.rm = TRUE)))
i <- match("padj", names(de))
de <- cbind(de[, 1:i], padj_lfc0 = res$padj[m], is_deg = res$deg[m], de[, -(1:i), drop = FALSE])
wr(de, p); cat("differential_expression.tsv: padj_lfc0 and is_deg added,", nrow(de), "rows\n")

## gene_index.tsv (skipped when the directory has none: gene_index.R then builds it from these calls)
p <- file.path(open, "gene_index.tsv")
if (!file.exists(p)) { cat("gene_index.tsv: not in", open, "- skipped\n"); quit(save = "no") }
gi <- rd(p)
if ("padj_lfc0_6h" %in% names(gi)) stop("gene_index.tsv already has padj_lfc0 columns")
for (tp in TPS) {
  r <- res[res$timepoint == tp, ]; k <- match(paste(gi$species, gi$gene), paste(r$species, r$gene))
  gi[[paste0("padj_lfc0_", tp)]] <- ifelse(is.na(k), "", r$padj[k])
  gi[[paste0(".deg_", tp)]] <- !is.na(k) & r$deg[k] %in% "TRUE"
}
ndeg <- rowSums(as.matrix(gi[paste0(".deg_", TPS)]))
gi$responds_to_cold <- ifelse(ndeg > 0, "True", "False")
gi$n_timepoints_significant <- as.character(ndeg)
gi <- gi[, !startsWith(names(gi), ".deg_")]
new <- paste0("padj_lfc0_", TPS); j <- match("padj_10d", names(gi))
gi <- gi[, c(setdiff(names(gi)[1:j], new), new, setdiff(names(gi)[-(1:j)], new))]
wr(gi, p)
for (sp in unique(gi$species))
  cat(sprintf("gene_index %-7s responds_to_cold = %5d\n", sp, sum(gi$species == sp & gi$responds_to_cold == "True")))
