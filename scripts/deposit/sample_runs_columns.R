#!/usr/bin/env Rscript
# sample_runs_columns.R — add each sample's ENA run accessions and sequencing batch to the deposit's samples.tsv.
#
# The values come from results/sample_ena_runs.tsv, written by scripts/deposit/sample_ena_runs.R (which traces every
# sample to its runs through the deposit's own sample-name and lane tables, and takes the platform from the read
# headers rather than from ENA's instrument_model field). Every sample must be matched exactly once.
#
# Edits in place, leaving every other column byte for byte as it was:
#   <open>/samples.tsv   adds ena_runs (";"-separated run accessions) and sequencing_batch
# Usage: Rscript scripts/deposit/sample_runs_columns.R <deposit data/open dir>

open <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(open), dir.exists(open))
rd <- function(p) read.delim(p, colClasses = "character", na.strings = character(0), check.names = FALSE, quote = "")
p <- file.path(open, "samples.tsv"); s <- rd(p)
if (any(c("ena_runs", "sequencing_batch") %in% names(s))) stop("samples.tsv already has the run columns")
r <- rd("results/sample_ena_runs.tsv")
k <- match(paste(s$species, s$sample_id), paste(r$species, r$sample_id))
if (anyNA(k) || anyDuplicated(k) || nrow(r) != nrow(s)) stop("samples.tsv and results/sample_ena_runs.tsv do not match one to one")
s$ena_runs <- r$ena_runs[k]; s$sequencing_batch <- r$batch[k]
write.table(s, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("samples.tsv: ena_runs and sequencing_batch added for %d samples (%d runs); batches: %s\n", nrow(s),
            sum(lengths(strsplit(s$ena_runs, ";"))),
            paste(sprintf("%s %d", names(table(s$sequencing_batch)), table(s$sequencing_batch)), collapse = ", ")))
