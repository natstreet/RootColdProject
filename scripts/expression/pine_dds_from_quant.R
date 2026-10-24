#!/usr/bin/env Rscript
# pine_dds_from_quant.R — the Scots pine DESeq2 object rebuilt from the deposit's per-run Salmon files with tximport
# (ledger R059 step 1, audit A6: the deposited object's avgTxLength assay differs from the quant files in zero-count
# cells). One run file per pine sample (data/quantification/lane_map.tsv); gene = transcript id without ".mRNA.<n>";
# tximport(type = "salmon") with its own length handling, then DESeqDataSetFromTximport with the deposited object's
# sample table, sample order and design. Genes: the rows of the deposited object.
# Checks: counts equal the deposited counts for every gene and sample (the deposited counts are each run's reads
# summed per gene and rounded); the script reports how many avgTxLength cells change, split by zero and non-zero count.
# Usage (working tree with the deposit's data/): Rscript scripts/expression/pine_dds_from_quant.R <out.rda>
suppressPackageStartupMessages({ library(DESeq2); library(tximport); library(data.table) })
out <- commandArgs(trailingOnly = TRUE)[1]
e <- new.env(); load("data/dds/dds_Ps.rda", envir = e); dep <- get(ls(e)[1], e)
lm <- fread("data/quantification/lane_map.tsv", colClasses = "character")[genotype == "pine"]
if (anyDuplicated(lm$deposited_sample_id)) stop("a pine sample has more than one run file; this script expects one")
files <- setNames(file.path("data/quantification/lanes", paste0(lm$quantification_file, ".quant.sf.gz")), lm$deposited_sample_id)
files <- files[colnames(dep)]
if (anyNA(files) || !all(file.exists(files))) stop("missing run files for: ", paste(colnames(dep)[is.na(files) | !file.exists(files)], collapse = " "))
tx <- fread(cmd = paste("gzip -dc", shQuote(files[1])), select = "Name")$Name
t2g <- data.frame(tx = tx, gene = sub("\\.mRNA\\.[0-9]+$", "", tx))
txi <- tximport(files, type = "salmon", tx2gene = t2g)
keep <- rownames(dep)
if (!all(keep %in% rownames(txi$counts))) stop("deposited genes missing from the quant files")
txi[c("abundance", "counts", "length")] <- lapply(txi[c("abundance", "counts", "length")], function(m) m[keep, , drop = FALSE])
dds <- DESeqDataSetFromTximport(txi, colData(dep), design(dep))
if (!identical(unname(counts(dds)), unname(counts(dep)))) stop("rebuilt counts differ from the deposited counts")
old <- assay(dep, "avgTxLength"); new <- assay(dds, "avgTxLength"); z <- counts(dep) == 0
ch <- abs(new - old) > 1e-6 * pmax(abs(old), 1)
cat(sprintf("pine: %d genes x %d samples; counts equal the deposited counts\n", nrow(dds), ncol(dds)))
cat(sprintf("avgTxLength cells changed: %d of %d (zero-count cells %d of %d; non-zero-count cells %d of %d)\n",
            sum(ch), length(ch), sum(ch & z), sum(z), sum(ch & !z), sum(!z)))
dds_Ps <- dds; save(dds_Ps, file = out)
cat("written:", out, "\n")
