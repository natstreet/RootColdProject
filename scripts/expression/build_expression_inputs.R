#!/usr/bin/env Rscript
# build_expression_inputs.R — the DESeq2 count objects and the co-expression input matrices of the six root
# cold-stress experiments, with every sequencing run of every sample included.
#
# Counts per species:
#   Col-0, aspen, birch, Norway spruce, Scots pine: the deposited DESeq2 objects, unchanged (they already sum every
#     run). Pine sample 112 has two sequencing runs, S2 and S9; its deposited column equals a quantification of a
#     merged read file (59.5 M reads), and which runs went into that file was not confirmed (ledger A033, N126).
#   Ost-0: the deposited object, with sample OstR1001 given both of its sequencing runs: each run's salmon read
#     counts summed per gene and rounded, then added over runs (the convention of the other samples, which the
#     script checks by reproducing their deposited counts).
# Co-expression input matrices: DESeq2 variance-stabilising transformation (blind = TRUE) minus the matrix's
# smallest value, as scripts/ComPlEx/Network_Generation_setup.R; columns named as in data/expression/sample_names.tsv.
#
# Run from the repository root:
#   Rscript scripts/expression/build_expression_inputs.R <deposit data/dds dir> <salmon dir with OstR1001 runs>
# Outputs: data/dds_rebuilt/<object>.rda (six DESeq2 objects) and data/expression/<species>_root_cold_stress_expression.txt
suppressPackageStartupMessages({ library(DESeq2); library(data.table) })
a <- commandArgs(trailingOnly = TRUE)
DDS <- a[1]; SALMON <- a[2]
dir.create("data/dds_rebuilt", showWarnings = FALSE)
load_dds <- function(f) { e <- new.env(); load(file.path(DDS, f), envir = e); get(ls(e)[1], e) }
names_tab <- read.delim("data/expression/sample_names.tsv", colClasses = "character")
names_tab$matrix_name[names_tab$matrix_name == "Ost.044"] <- "Ost4.04"     # typo in the original Ost-0 matrix

lane_genes <- function(dir) {                   # one run's read counts summed per gene, rounded
  q <- fread(file.path(dir, "quant.sf"), select = c("Name", "NumReads"))
  q[, gene := sub("\\.[0-9]+$", "", Name)]
  g <- q[, .(n = round(sum(NumReads))), by = gene]; setNames(g$n, g$gene)
}
runs_of <- function(sample) list.dirs(SALMON, recursive = FALSE)[startsWith(basename(list.dirs(SALMON, recursive = FALSE)), paste0(sample, "_"))]
sample_counts <- function(sample, genes) {
  m <- vapply(runs_of(sample), function(r) { x <- lane_genes(r)[genes]; x[is.na(x)] <- 0; x }, numeric(length(genes)))
  rowSums(matrix(m, nrow = length(genes)))
}

objects <- list()
for (s in c("col0", "aspen", "birch", "spruce", "pine")) {
  f <- c(col0 = "dds_col0_soil.rda", aspen = "dds_Pt_new.rda", birch = "dds_Bp_new.rda", spruce = "dds_Pa.rda", pine = "dds_Ps.rda")[[s]]
  objects[[s]] <- list(file = f, dds = load_dds(f))
}

# Ost-0: OstR1001 with both runs
d <- load_dds("dds_ost0_soil.rda"); m <- counts(d)
for (chk in c("OstRWT01", "OstR0601", "OstR2401")) {    # the run convention reproduces deposited single-run samples
  if (!identical(unname(sample_counts(chk, rownames(m))), unname(as.numeric(m[, chk])))) stop("run convention does not reproduce ", chk)
}
runs <- runs_of("OstR1001"); if (length(runs) != 2) stop("OstR1001: expected two runs, found ", length(runs))
new1001 <- sample_counts("OstR1001", rownames(m))
cat(sprintf("Ost-0 OstR1001: %d runs; deposited %.0f reads, both runs %.0f reads\n", length(runs), sum(m[, "OstR1001"]), sum(new1001)))
m[, "OstR1001"] <- new1001; storage.mode(m) <- "integer"
d2 <- DESeqDataSetFromMatrix(m, colData(d), design(d))
objects$ost0 <- list(file = "dds_ost0_soil.rda", dds = d2)

species_file <- c(col0 = "col0", ost0 = "ost0", aspen = "aspen", birch = "birch", spruce = "spruce", pine = "pine")
for (s in names(objects)) {
  dd <- objects[[s]]$dds
  nm <- sub("\\.rda$", "", objects[[s]]$file); assign(nm, dd)
  save(list = nm, file = file.path("data/dds_rebuilt", objects[[s]]$file))
  v <- assay(varianceStabilizingTransformation(dd, blind = TRUE)); v <- v - min(v)
  t <- names_tab[names_tab$species == s, ]
  if (!setequal(t$dds_sample, colnames(v))) stop(s, ": samples differ from sample_names.tsv")
  v <- v[, t$dds_sample]; colnames(v) <- t$matrix_name
  write.table(data.frame(Genes = rownames(v), v, check.names = FALSE),
              sprintf("data/expression/%s_root_cold_stress_expression.txt", species_file[[s]]),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%-6s %d genes x %d samples -> data/expression/%s_root_cold_stress_expression.txt\n", s, nrow(v), ncol(v), species_file[[s]]))
}
