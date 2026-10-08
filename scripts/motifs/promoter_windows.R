#!/usr/bin/env Rscript
# promoter_windows.R — promoter windows anchored on the start codon, for the promoter motif track (DECISIONS.md D062
# item 4, D063; A035: "1 kb upstream of the ATG (2 kb sensitivity), not the annotated TSS").
#
# For each listed gene: the transcript with the longest CDS (ties: first transcript ID in sort order); its start codon is
# the first CDS base on the gene's strand; the window is the WIN bases immediately upstream of it, on that strand, clipped
# at sequence ends. Also recorded per window (A035 covariates): bases overlapping any other gene's span (either strand),
# and, when a repeat interval file is given, bases overlapping repeats (union of intervals).
# Inputs: <gff3> (gene / mRNA / CDS features; gz allowed), <genes.txt>, <fasta .fai>, <seq prefix> (prepended to GFF
#   seqids to give FASTA names, "-" for none), <repeat intervals> (BED or GFF, gz allowed, "-" for none), <out prefix>
# Outputs: <out prefix>_<WIN>.bed (BED6, name = gene ID) and <out prefix>_windows.tsv (one row per gene and window)
# Stops if a listed gene has no CDS or its sequence is not in the .fai.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
GFF <- a[1]; GENES <- a[2]; FAI <- a[3]; PFX <- if (a[4] == "-") "" else a[4]; REP <- a[5]; OUT <- a[6]
WINS <- c(1000L, 2000L)
nocomment <- function(f) paste(if (grepl("\\.gz$", f)) "gzip -dc" else "cat", shQuote(f), "| grep -v '^#'")   # comment lines confuse fread
gff <- fread(cmd = nocomment(GFF), sep = "\t", header = FALSE, quote = "", fill = TRUE, select = 1:9, colClasses = "character")
setnames(gff, c("seq", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
gff[, `:=`(start = as.integer(start), end = as.integer(end))]
att <- function(x, k) sub(sprintf(".*(^|;)%s=([^;]*).*", k), "\\2", ifelse(grepl(sprintf("(^|;)%s=", k), x), x, NA))
genes <- gff[type == "gene", .(gene = att(attr, "ID"), seq, gstart = start, gend = end, strand)]
mrna <- gff[type %in% c("mRNA", "transcript"), .(tx = att(attr, "ID"), gene = att(attr, "Parent"))]
cds <- gff[type == "CDS", .(parent = att(attr, "Parent"), start, end)]
cds <- cds[, .(tx = unlist(strsplit(parent, ","))), by = .(start, end, parent)]   # Araport lists several parents
cds <- merge(cds, mrna, by = "tx")
want <- readLines(GENES)
txl <- cds[gene %in% want, .(cds_len = sum(end - start + 1L), cmin = min(start), cmax = max(end)), by = .(gene, tx)]
setorder(txl, gene, -cds_len, tx); txl <- txl[, .SD[1], by = gene]
miss <- setdiff(want, txl$gene); if (length(miss)) stop(length(miss), " listed genes have no CDS, e.g. ", paste(head(miss, 5), collapse = " "))
w0 <- merge(txl, genes, by = "gene")
fai <- fread(FAI, header = FALSE, select = 1:2, col.names = c("fseq", "len"))
w0[, fseq := paste0(PFX, seq)]
if (w0[!fseq %in% fai$fseq, .N]) stop("sequence names not in the .fai: ", paste(head(unique(w0[!fseq %in% fai$fseq, fseq])), collapse = " "))
w0 <- merge(w0, fai, by = "fseq")
w0[, atg := ifelse(strand == "+", cmin, cmax)]
rpt <- NULL
if (REP != "-") {
  isg <- grepl("\\.gff3?(\\.gz)?$", REP)
  r <- fread(cmd = nocomment(REP), sep = "\t", header = FALSE, quote = "", fill = TRUE, select = if (isg) c(1L, 4L, 5L) else 1:3,
             colClasses = "character", col.names = c("seq", "a", "b"))
  rpt <- if (isg) r[, .(seq, s = as.integer(a), e = as.integer(b))] else r[, .(seq, s = as.integer(a) + 1L, e = as.integer(b))]
  rpt <- rpt[!is.na(s) & seq %in% unique(w0$seq)]; setkey(rpt, seq, s, e); cat("repeat intervals on gene sequences:", nrow(rpt), "\n")
}
union_bp <- function(n, x) {   # per window id 1..n: bases covered by the union of its intervals x (id, a, b)
  out <- integer(n); if (!nrow(x)) return(out)
  u <- x[order(id, a)][, { m <- cummax(b); g <- cumsum(c(TRUE, a[-1] > m[-.N] + 1L)); .(bp = sum(tapply(b, g, max) - tapply(a, g, min) + 1L)) }, by = id]
  out[u$id] <- as.integer(u$bp); out }
res <- rbindlist(lapply(WINS, function(W) {
  w <- copy(w0)[, `:=`(win = W, id = .I)]
  w[strand == "+", `:=`(ws = pmax(1L, atg - W), we = atg - 1L)]
  w[strand == "-", `:=`(ws = atg + 1L, we = pmin(len, atg + W))]
  w[, win_len := pmax(0L, we - ws + 1L)]
  og <- genes[, .(seq, s = gstart, e = gend, other = gene)]
  setkey(og, seq, s, e)
  ov <- foverlaps(w[win_len > 0, .(id, gene, seq, s = ws, e = we)], og, by.x = c("seq", "s", "e"), nomatch = NULL)[other != gene]
  w[, neighbour_gene_bp := union_bp(nrow(w), ov[, .(id, a = pmax(s, i.s), b = pmin(e, i.e))])]
  if (!is.null(rpt)) {
    rv <- foverlaps(w[win_len > 0, .(id, seq, s = ws, e = we)], rpt, by.x = c("seq", "s", "e"), nomatch = NULL)
    w[, repeat_bp := union_bp(nrow(w), rv[, .(id, a = pmax(s, i.s), b = pmin(e, i.e))])]
  } else w[, repeat_bp := NA_integer_]
  fwrite(w[win_len > 0, .(fseq, ws - 1L, we, gene, 0L, strand)], sprintf("%s_%d.bed", OUT, W), sep = "\t", col.names = FALSE)
  w[, .(gene, win = W, transcript = tx, seq = fseq, strand, atg, ws, we, win_len, clipped = win_len < W, neighbour_gene_bp, repeat_bp)] }))
fwrite(res, sprintf("%s_windows.tsv", OUT), sep = "\t")
print(res[, .(genes = .N, clipped = sum(clipped), mean_len = round(mean(win_len)), mean_neighbour_bp = round(mean(neighbour_gene_bp)),
              mean_repeat_bp = if (all(is.na(repeat_bp))) NA_real_ else round(mean(repeat_bp))), by = win], row.names = FALSE)
