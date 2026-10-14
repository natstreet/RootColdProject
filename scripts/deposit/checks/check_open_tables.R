#!/usr/bin/env Rscript
# check_open_tables.R — the open tables of data/open/ against the objects they are exported from, and the sample table
# against the deposited DESeq2 objects. Writes one row per check with ok = TRUE/FALSE; stops when any check fails.
#   co_expressologs.tsv    = data/ComPlEx/co_expressologs.RData (same pairs, MaxpVal and effectSize within 1e-12 relative)
#   pairwise_comparisons.tsv = the 15 data/ComPlEx/PairwiseComparisons/*.RData tables (rows and p-values; the RData tables
#                            store the p-value columns as text)
#   conserved_core.tsv     = the orthogroups of data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS
#   samples.tsv            every sample of every DESeq2 object, once, through data/quantification/sample_name_map.tsv, with
#                          the DESeq2 object's time point (its design variable, "_5C" dropped)
# Usage: Rscript check_open_tables.R <deposit root> <out.tsv>
suppressMessages({ library(DESeq2); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); R <- a[1]; OUT <- a[2]
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], e) }
rows <- list(); add <- function(check, ok, detail) rows[[length(rows) + 1]] <<- data.table(check = check, ok = ok, detail = detail)
rel <- function(x, y) max(abs(x - y) / pmax(abs(y), 1e-300), na.rm = TRUE)

co <- as.data.table(ld(file.path(R, "data/ComPlEx/co_expressologs.RData"))); ct <- fread(file.path(R, "data/open/co_expressologs.tsv"))
k1 <- paste(co$Species1, co$GeneSpecies1, co$Species2, co$GeneSpecies2, co$OrthologGroup)
k2 <- paste(ct$Species1, ct$GeneSpecies1, ct$Species2, ct$GeneSpecies2, ct$OrthologGroup)
ok <- nrow(co) == nrow(ct) && setequal(k1, k2)
d <- if (ok) { m <- match(k1, k2); max(rel(ct$MaxpVal[m], as.numeric(co$MaxpVal)), rel(ct$effectSize[m], co$effectSize)) } else NA
add("co_expressologs.tsv = RData", ok && d < 1e-12, sprintf("rows %d / %d; max relative difference %s", nrow(ct), nrow(co), format(d)))

pw <- fread(file.path(R, "data/open/pairwise_comparisons.tsv")); n_obj <- 0; worst <- 0; keys_ok <- TRUE
for (f in list.files(file.path(R, "data/ComPlEx/PairwiseComparisons"), full.names = TRUE)) {
  x <- as.data.table(ld(f)); sp <- strsplit(sub("^comparison-(.*)-abspearson.*$", "\\1", basename(f)), "-")[[1]]; n_obj <- n_obj + nrow(x)
  y <- pw[species_pair == paste(sp, collapse = "-")]
  kx <- paste(x$Species1, x$Species2, x$OrthoGroup); ky <- paste(y$Species1, y$Species2, y$OrthoGroup)
  if (nrow(x) != nrow(y) || !setequal(kx, ky)) { keys_ok <- FALSE; next }
  m <- match(kx, ky)
  for (cc in c("Species1.p.val", "Species2.p.val", "Max.p.val")) worst <- max(worst, rel(y[[cc]][m], as.numeric(x[[cc]])))
}
add("pairwise_comparisons.tsv = 15 RData tables", keys_ok && nrow(pw) == n_obj && worst < 1e-12,
    sprintf("rows %d / %d; keys equal %s; max relative p-value difference %s", nrow(pw), n_obj, keys_ok, format(worst)))

cc <- fread(file.path(R, "data/open/conserved_core.tsv")); core <- as.data.frame(readRDS(file.path(R, "data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")))
add("conserved_core.tsv = core RDS", nrow(cc) == nrow(core) && setequal(cc$OrthoGroup, core$OrthoGroup), sprintf("rows %d / %d", nrow(cc), nrow(core)))

G <- c(col0 = "dds_col0_soil.rda", ost0 = "dds_ost0_soil.rda", aspen = "dds_Pt_new.rda", birch = "dds_Bp_new.rda", spruce = "dds_Pa.rda", pine = "dds_Ps.rda")
sm <- fread(file.path(R, "data/open/samples.tsv"), colClasses = "character"); nm <- fread(file.path(R, "data/quantification/sample_name_map.tsv"), colClasses = "character")
for (g in names(G)) {
  dds <- ld(file.path(R, "data/dds", G[[g]])); v <- all.vars(design(dds))[1]
  lvl <- sub("_5C$", "", as.character(colData(dds)[[v]])); pub <- setNames(nm[genotype == g]$deposited_sample_id, nm[genotype == g]$sequencing_name)[colnames(dds)]
  s <- sm[species == g]; tp <- setNames(s$timepoint, s$sample_id)[pub]
  ok <- !anyNA(pub) && nrow(s) == ncol(dds) && setequal(s$sample_id, pub) && all(tp == lvl)
  add(sprintf("samples.tsv = DESeq2 samples and time points (%s)", g), ok, sprintf("samples %d / %d; time points differing %d", nrow(s), ncol(dds), sum(tp != lvl, na.rm = TRUE)))
}
out <- rbindlist(rows); fwrite(out, OUT, sep = "\t"); print(out)
if (!all(out$ok)) stop("open-table checks failed")
