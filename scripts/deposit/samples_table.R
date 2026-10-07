#!/usr/bin/env Rscript
# samples_table.R — the deposit's samples.tsv (design columns): one row per column of the six expression matrices
# (data/expression/<species>_root_cold_stress_expression.txt), decoded from the column name <prefix><replicate>.<code>:
# code 01 = control, 02 = 6h, 03 = 24h, 04 = 3d, 05 = 10d (5 °C for every cold timepoint). Rows are ordered by species,
# then timepoint as text, then replicate. scripts/deposit/sample_runs_columns.R then adds ena_runs and sequencing_batch.
# Usage: Rscript scripts/deposit/samples_table.R <deposit data dir> <out samples.tsv>
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); D <- a[1]
SP <- data.table(species = c("col0", "ost0", "aspen", "birch", "spruce", "pine"),
                 species_binomial = c("Arabidopsis thaliana", "Arabidopsis thaliana", "Populus tremula", "Betula pendula",
                                      "Picea abies", "Pinus sylvestris"),
                 ecotype = c("Col-0", "Ost-0", "", "", "", ""), prefix = c("Col", "Ost", "Pt", "Bp", "Pa", "Ps"))
TP <- data.table(code = sprintf("%02d", 1:5), timepoint = c("control", "6h", "24h", "3d", "10d"), hours_cold = c(0L, 6L, 24L, 72L, 240L))
s <- rbindlist(lapply(seq_len(nrow(SP)), function(i) {
  id <- setdiff(names(fread(file.path(D, "expression", sprintf("%s_root_cold_stress_expression.txt", SP$species[i])), nrows = 0)), "Genes")
  m <- regmatches(id, regexec(sprintf("^%s([0-9]+)[.]([0-9]{2})$", SP$prefix[i]), id))
  if (any(lengths(m) != 3)) stop(SP$species[i], ": column names not of the form <prefix><replicate>.<code>")
  data.table(species = SP$species[i], species_binomial = SP$species_binomial[i], ecotype = SP$ecotype[i], sample_id = id,
             replicate = as.integer(vapply(m, `[`, "", 2)), code = vapply(m, `[`, "", 3))
}))
s <- merge(s, TP, by = "code", sort = FALSE)
if (anyNA(s$timepoint)) stop("unknown timepoint code")
s[, condition := ifelse(timepoint == "control", "control", "cold")][, temperature_C := ifelse(timepoint == "control", "", "5")]
setorder(s, species, timepoint, replicate)
s <- s[, .(species, species_binomial, ecotype, sample_id, replicate, timepoint, hours_cold, condition, temperature_C)]
fwrite(s, a[2], sep = "\t", quote = FALSE)
cat("samples.tsv:", nrow(s), "samples;", paste(s[, .N, by = species][, sprintf("%s %d", species, N)], collapse = ", "), "\n")
