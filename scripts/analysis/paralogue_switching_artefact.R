#!/usr/bin/env Rscript
# paralogue_switching_artefact.R — the quantification-artefact check of each candidate genuine switch of
# paralogue_switching.R (A035, pre-specified): Ost-0 reads were quantified against the Col-0 (Araport11) index, so
# accession polymorphisms can move multi-mapping reads between near-identical paralogues. Per candidate: CDS identity of
# the two loci (EMBOSS needle global alignment of their .1 CDS, gap open 10, extend 0.5; identity over the alignment
# length), whether they are tandem neighbours (only the different-chromosome case is decided here, as not tandem; a
# same-chromosome pair stops the script, since that needs the full gene order),
# and their expression correlation within Col-0 and within Ost-0. CDS identity >= 95% -> "possible quantification
# artefact", otherwise the candidate stays a candidate genuine switch.
# Inputs (working tree): results/paralogue_switching.tsv, data/paralogue_switching/{candidate_cds.fa,candidate_positions.tsv}
#   (extracted from Araport11, SOURCE.txt), data/expression/{col0,ost0}_root_cold_stress_expression.txt; needle on PATH or via pixi
# Output: results/paralogue_switching_artefact.tsv. Usage (working tree): Rscript scripts/analysis/paralogue_switching_artefact.R
suppressPackageStartupMessages(library(data.table))
x <- fread("results/paralogue_switching.tsv")[category == "candidate genuine switch"]
fa <- readLines("data/paralogue_switching/candidate_cds.fa"); pos <- fread("data/paralogue_switching/candidate_positions.tsv", header = FALSE,
  col.names = c("gene", "seq", "start", "end", "strand"))
cds1 <- function(g) { h <- grep(paste0("^>", g, "\\.1 "), fa); if (!length(h)) stop("no .1 CDS for ", g)
  e <- c(grep("^>", fa), length(fa) + 1); e <- e[e > h][1] - 1; c(paste0(">", g, ".1"), fa[(h + 1):e]) }
needle <- if (nzchar(Sys.which("needle"))) "needle" else "pixi exec --spec emboss -c bioconda -c conda-forge -- needle"
mexp <- function(g) { y <- fread(sprintf("data/expression/%s_root_cold_stress_expression.txt", g)); m <- as.matrix(y[, -1]); rownames(m) <- y[[1]]; m }
mc <- mexp("col0"); mo <- mexp("ost0")
res <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
  a <- x$col0_member[i]; b <- x$ost0_member[i]; td <- tempfile(); dir.create(td)
  writeLines(cds1(a), file.path(td, "a.fa")); writeLines(cds1(b), file.path(td, "b.fa"))
  system(sprintf("%s -asequence %s -bsequence %s -gapopen 10 -gapextend 0.5 -outfile %s -auto", needle, file.path(td, "a.fa"),
                 file.path(td, "b.fa"), file.path(td, "n.txt")), ignore.stdout = TRUE, ignore.stderr = TRUE)
  idl <- grep("^# Identity:", readLines(file.path(td, "n.txt")), value = TRUE)
  ident <- as.numeric(sub(".*\\(\\s*([0-9.]+)%\\).*", "\\1", idl))
  pa <- pos[gene == a]; pb <- pos[gene == b]
  if (!nrow(pa) || !nrow(pb)) stop("no position for ", a, " or ", b)
  if (pa$seq == pb$seq) stop(a, " and ", b, " are on the same sequence: the tandem check needs the full gene order")
  tandem <- FALSE
  data.table(OrthoGroup = x$OrthoGroup[i], col0_member = a, ost0_member = b, cds_identity_pct = ident, alignment = sub("^# Identity:\\s*", "", idl),
             same_sequence = FALSE, tandem_neighbours = tandem,
             r_within_col0 = round(cor(mc[a, ], mc[b, ]), 3), r_within_ost0 = round(cor(mo[a, ], mo[b, ]), 3),
             call = ifelse(ident >= 95, "possible quantification artefact", "candidate genuine switch")) }))
fwrite(res, "results/paralogue_switching_artefact.tsv", sep = "\t", na = ""); print(res)
