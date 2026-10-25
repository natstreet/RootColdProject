#!/usr/bin/env Rscript
# hr_family_domains.R — protein-domain identity of the families behind "plant-type hypersensitive response" (GO:0009626) among the
# conifer-only co-expressologs (ledger R066, D091).
# Families: the union of the GO:0009626 family_ids of the R064 family-level tables (Norway spruce, Scots pine). Proteins: every Norway
# spruce (Picea_abies) and Scots pine (Pinus_sylvestris) gene of these orthogroups (Orthogroups_20240613.tsv), one protein per gene
# from the reference proteomes' longest-isoform files (header attribute gene=). hmmscan against Pfam-A with --cut_ga; a protein
# carries a domain when a row of the per-sequence table (--tblout) for that Pfam accession has an inclusion count (column 18, inc)
# of at least 1, read as an integer.
# Domain sets: TIR PF01582, PF13676; NB-ARC PF00931; LRR PF00560, PF07725, PF13855, PF12799; RPW8 PF05659; kinase PF00069, PF07714.
# Gene class (first rule that applies): RNL = NB-ARC and RPW8; TNL = NB-ARC and TIR; CNL-like = NB-ARC and LRR; NB-only = NB-ARC;
#   TIR-only = TIR; RLK = kinase and LRR; other. Family class: the most frequent gene class over its spruce and pine genes ("other"
#   included; ties named together). NLR-class family: family class TNL, CNL-like, RNL or NB-only.
# Sensitivity (class_ext, nlr_class_ext): Pfam 37+ models NLR regions next to NB-ARC as separate families, which the sets above
#   predate. NB evidence extended with the disease-resistance winged-helix domains WHD_DRP PF23559 and WHD_ROQ1 PF23282; LRR set
#   extended with LRR_14 PF23598, LRR_R13L1-DRL21 PF25019, LRR_RPS2 PF23247 and LRR_6 PF13516. Same class rules.
# The description line per family is assembled from the data: family class, gene-class counts, the most frequent Pfam domains
#   (genes carrying them), the number of the 27 OrthoFinder species with genes, and the Arabidopsis members' symbols
#   (gene_aliases_20140331.txt, first symbol per locus).
# Also: the GO:0009626 annotated study genes of each conifer (R064 gene_level_study_genes_<conifer>.tsv) in NLR-class families.
# Inputs: <R064 dir> (conifer_only_go.R output), data/ComPlEx/Orthogroups_20240613.tsv, data/annotation/gene_aliases_20140331.txt,
#   <pressed Pfam-A.hmm>, <Picab02 longest aa .fa.gz>, <Pinsy01 longest aa .fa.gz>
# Outputs (out dir): proteins.fa, hmmscan_tblout.txt, hmmscan_domtblout.txt, gene_domains.tsv, family_domains.tsv,
#   hr_genes_in_nlr_families.tsv, summary.txt
# Usage (working tree): Rscript scripts/analysis/hr_family_domains.R <out dir> <R064 dir> <Pfam-A.hmm> <spruce aa> <pine aa> [cpu]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); OUT <- a[1]; R064 <- a[2]; PFAM <- a[3]; AA <- c(spruce = a[4], pine = a[5])
CPU <- if (length(a) >= 6) a[6] else "4"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TERM <- "GO:0009626"
DOM <- list(TIR = c("PF01582", "PF13676"), NB_ARC = "PF00931", LRR = c("PF00560", "PF07725", "PF13855", "PF12799"),
            RPW8 = "PF05659", kinase = c("PF00069", "PF07714"))
DOM_EXT <- modifyList(DOM, list(NB_ARC = c(DOM$NB_ARC, "PF23559", "PF23282"), LRR = c(DOM$LRR, "PF23598", "PF25019", "PF23247", "PF13516")))
classify <- function(NB_ARC, TIR, LRR, RPW8, kinase)
  fifelse(NB_ARC & RPW8, "RNL", fifelse(NB_ARC & TIR, "TNL", fifelse(NB_ARC & LRR, "CNL-like", fifelse(NB_ARC, "NB-only",
  fifelse(TIR, "TIR-only", fifelse(kinase & LRR, "RLK", "other"))))))
# known-answer test: RPW8 outranks TIR (an NB-ARC+TIR+RPW8 protein is RNL), TIR without NB-ARC is TIR-only, kinase without LRR is other
stopifnot(identical(classify(c(T, T, T, T, F, F, F, F), c(T, T, F, F, T, F, F, F), c(F, T, T, F, F, T, F, T), c(T, F, F, F, F, F, F, F),
                             c(F, F, F, F, F, T, T, F)),
                    c("RNL", "TNL", "CNL-like", "NB-only", "TIR-only", "RLK", "other", "other")))
COL <- c(spruce = "Picea_abies", pine = "Pinus_sylvestris")
NLR <- c("TNL", "CNL-like", "RNL", "NB-only")

## families and HR study genes (R064)
fam <- sort(unique(unlist(lapply(names(COL), function(cf) {
  f <- fread(file.path(R064, sprintf("family_level_%s.tsv", cf)))[GO.ID == TERM]
  strsplit(f$family_ids, ",", fixed = TRUE)[[1]] }))))
hr <- rbindlist(lapply(names(COL), function(cf) fread(file.path(R064, sprintf("gene_level_study_genes_%s.tsv", cf)))[GO.ID == TERM][, conifer := cf]))

## genes of each family
og <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE, colClasses = "character")   # header one field short:
names(og)[1] <- "OrthoGroup"; og <- as.data.table(og)[OrthoGroup %in% fam]; stopifnot(nrow(og) == length(fam), ncol(og) == 28)  # col 1 = row names
spl <- function(x) { x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]]); x[nzchar(x)] }
genes <- rbindlist(lapply(names(COL), function(cf) og[, .(gene = spl(get(COL[[cf]]))), by = OrthoGroup][, conifer := cf]))
miss <- hr[!genes, on = .(gene, OrthoGroup)]
if (nrow(miss)) stop("HR study genes not among their orthogroup's genes: ", paste(miss$gene, collapse = ", "))

## proteins: one per gene, from the longest-isoform proteomes
read_fa <- function(f) {
  l <- readLines(gzfile(f)); h <- grep("^>", l)
  id <- sub("^>(\\S+).*$", "\\1", l[h]); g <- sub("^.*\\bgene=(\\S+).*$", "\\1", l[h])
  sq <- vapply(seq_along(h), function(i) paste(l[seq.int(h[i] + 1, c(h[-1] - 1, length(l))[i])], collapse = ""), "")
  data.table(protein = id, gene = g, seq = sq)
}
pr <- rbindlist(lapply(names(AA), function(cf) read_fa(AA[[cf]])[gene %in% genes[conifer == cf, gene]]))
stopifnot(!anyDuplicated(pr$gene))
nf <- setdiff(genes$gene, pr$gene); if (length(nf)) stop("genes without a protein: ", paste(head(nf, 10), collapse = ", "))
fa <- file.path(OUT, "proteins.fa"); writeLines(paste0(">", pr$protein, "\n", pr$seq), fa)

## hmmscan
tbl <- file.path(OUT, "hmmscan_tblout.txt"); dom <- file.path(OUT, "hmmscan_domtblout.txt")
st <- system2("hmmscan", c("--cut_ga", "--cpu", CPU, "--tblout", tbl, "--domtblout", dom, "-o", "/dev/null", PFAM, fa))
if (st != 0) stop("hmmscan failed: ", st)
x <- readLines(tbl); x <- x[!startsWith(x, "#")]
f <- strsplit(trimws(x), "\\s+")
hits <- data.table(pfam = vapply(f, `[`, "", 1), acc = sub("\\..*$", "", vapply(f, `[`, "", 2)), protein = vapply(f, `[`, "", 3),
                   rep = as.integer(vapply(f, `[`, "", 17)), inc = as.integer(vapply(f, `[`, "", 18)))
stopifnot(!anyNA(hits$inc), all(hits$inc <= hits$rep), all(hits$protein %in% pr$protein))   # inc: included domains, <= reported
hits <- hits[inc >= 1][, gene := pr$gene[match(protein, pr$protein)]]

## gene and family classes
gd <- merge(genes, pr[, .(gene, protein)], by = "gene")
for (d in names(DOM)) gd[, (d) := gene %in% hits[acc %in% DOM[[d]], gene]]
gd[, class := classify(NB_ARC, TIR, LRR, RPW8, kinase)]
ext <- lapply(DOM_EXT, function(v) gd$gene %in% hits[acc %in% v, gene])
gd[, class_ext := classify(ext$NB_ARC, ext$TIR, ext$LRR, ext$RPW8, ext$kinase)]
fwrite(gd[order(OrthoGroup, conifer, gene)], file.path(OUT, "gene_domains.tsv"), sep = "\t")
al <- fread("data/annotation/gene_aliases_20140331.txt", quote = "")[, .SD[1], by = locus_name]
ath <- og[, .(at = spl(Arabidopsis_thaliana)), by = OrthoGroup][, locus := toupper(sub("\\.[0-9]+$", "", at))]
ath[, symbol := al$symbol[match(locus, al$locus_name)]]
nsp <- og[, .(n_species = sum(vapply(.SD, function(v) length(spl(v)) > 0, TRUE))), by = OrthoGroup, .SDcols = 2:ncol(og)]
fd <- rbindlist(lapply(fam, function(o) {
  g <- gd[OrthoGroup == o]; row <- list(OrthoGroup = o)
  for (cf in names(COL)) { h <- g[conifer == cf]; row[[paste0(cf, "_genes")]] <- nrow(h)
    for (d in names(DOM)) row[[paste0(cf, "_", d)]] <- sum(h[[d]]) }
  tc <- sort(table(g$class), decreasing = TRUE); top <- names(tc)[tc == max(tc)]
  row$class <- paste(top, collapse = "/"); row$nlr_class <- all(top %in% NLR)
  row$class_counts <- paste(sprintf("%s %d", names(tc), tc), collapse = ", ")
  te <- sort(table(g$class_ext), decreasing = TRUE); tp <- names(te)[te == max(te)]
  row$class_ext <- paste(tp, collapse = "/"); row$nlr_class_ext <- all(tp %in% NLR)
  row$class_ext_counts <- paste(sprintf("%s %d", names(te), te), collapse = ", ")
  pf <- hits[gene %in% g$gene, .(genes = uniqueN(gene)), by = pfam][order(-genes, pfam)]
  sy <- ath[OrthoGroup == o]; sy <- unique(fifelse(is.na(sy$symbol) | sy$symbol == sy$locus, sy$locus, paste0(sy$symbol, " (", sy$locus, ")")))
  row$description <- sprintf("%s; genes by class: %s; Pfam (genes): %s; species with genes: %d of 27; Arabidopsis: %s",
    row$class, row$class_counts, if (nrow(pf)) paste(sprintf("%s (%d)", head(pf$pfam, 5), head(pf$genes, 5)), collapse = ", ") else "none",
    nsp[OrthoGroup == o, n_species], if (length(sy)) paste(c(head(sy, 6), if (length(sy) > 6) sprintf("and %d more", length(sy) - 6)), collapse = ", ") else "none")
  as.data.table(row) }))
fwrite(fd, file.path(OUT, "family_domains.tsv"), sep = "\t")

## HR study genes in NLR-class families
hr[, `:=`(nlr_family = OrthoGroup %in% fd[nlr_class == TRUE, OrthoGroup], nlr_family_ext = OrthoGroup %in% fd[nlr_class_ext == TRUE, OrthoGroup])]
hr <- merge(hr, gd[, .(gene, gene_class = class, gene_class_ext = class_ext)], by = "gene", all.x = TRUE)
fwrite(hr[order(conifer, OrthoGroup, gene)], file.path(OUT, "hr_genes_in_nlr_families.tsv"), sep = "\t")
lines <- c(sprintf("families: %d; proteins scanned: %d (spruce %d, pine %d); Pfam-A models: %s", length(fam), nrow(pr),
                   gd[conifer == "spruce", .N], gd[conifer == "pine", .N], system2("grep", c("-c", "'^NAME'", PFAM), stdout = TRUE)),
           sprintf("NLR-class families: %d of %d; family classes: %s", fd[nlr_class == TRUE, .N], nrow(fd),
                   paste(sprintf("%s %d", names(table(fd$class)), table(fd$class)), collapse = ", ")),
           sprintf("%s annotated study genes in NLR-class families: %s", TERM,
                   paste(hr[, sprintf("%s %d of %d", conifer[1], sum(nlr_family), .N), by = conifer]$V1, collapse = "; ")),
           sprintf("sensitivity (extended NB and LRR sets): NLR-class families %d of %d; family classes: %s; study genes in NLR-class families: %s",
                   fd[nlr_class_ext == TRUE, .N], nrow(fd), paste(sprintf("%s %d", names(table(fd$class_ext)), table(fd$class_ext)), collapse = ", "),
                   paste(hr[, sprintf("%s %d of %d", conifer[1], sum(nlr_family_ext), .N), by = conifer]$V1, collapse = "; ")),
           sprintf("genes with a protein kinase domain: %d", gd[kinase == TRUE, .N]))
writeLines(lines, file.path(OUT, "summary.txt")); cat(lines, sep = "\n"); print(fd[, .(OrthoGroup, spruce_genes, pine_genes, class, class_ext)])
