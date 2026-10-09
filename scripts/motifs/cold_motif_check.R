#!/usr/bin/env Rscript
# cold_motif_check.R — are CBF binding sites (DRE/CRT) and the other cold-panel motifs over-represented in the promoters of
# the conserved core, against each genotype's expressed genes? (NRS 2026-10-09: firm up whether the core is CBF-driven in the
# conifers; the cold panel of A036.) Descriptive enrichment, criteria fixed before the run:
#   promoters    1 kb upstream of the ATG (promoter_windows.R, extract_promoters.py on micro); genes without an annotated start
#                codon excluded (D064); both strands scanned for the IUPAC consensus (cold_panel_motifs.tsv)
#   foreground   the conserved core's representative member in that genotype (133 orthogroups, D066)
#   background   the other expressed genes of the genotype (ComPlEx filters) with a promoter
#   test         per genotype x motif: Fisher's exact test on presence (at least one site); odds ratio with 95% CI; BH over
#                genotypes x motifs. Enriched if q < 0.05 and OR > 1. A null is read through the CI's upper bound.
#   control      Col-0: the CBF-activated genes of Zhao et al. (2016, Table S3) must be enriched for DRE/CRT against the same
#                background (q < 0.05, OR > 1), otherwise the scan is broken and the script stops; also our per-gene DRE/CRT
#                presence is compared with Table S3's "Number of CRT/DRE motif" column (agreement reported).
#   CBF genes    descriptive: the loci with symbol CBF1, CBF2 or CBF3 in data/annotation/gene_aliases_20140331.txt, their
#                orthogroup and its core_status (orthogroup_core_status from core_status_columns.R)
# Inputs: <promoter dir> (<genotype>_1000.fa, start_codons.tsv), cold_panel_motifs.tsv, <external sets tsv>,
#   data/external_gene_sets/PP2016-00533R1_Table_S3.xls, results/core_stability_label.tsv, data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS,
#   <orthogroup_core_status.tsv>, data/ComPlEx/Orthogroups_20240613.tsv (via the comparison tables), data/annotation/gene_aliases_20140331.txt
# Outputs (<out dir>): cold_motif_enrichment.tsv, cold_motif_control.tsv, cbf_genes_status.tsv
# Usage (working tree): Rscript cold_motif_check.R <promoter dir> <motifs tsv> <external sets tsv> <orthogroup_core_status.tsv> <out dir>
suppressPackageStartupMessages({ library(data.table); library(readxl) })
a <- commandArgs(trailingOnly = TRUE); PD <- a[1]; MT <- a[2]; ES <- a[3]; CS <- a[4]; OUT <- a[5]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
G <- c("col0", "aspen", "birch", "spruce", "pine")
mot <- fread(MT)
iu <- c(A = "A", C = "C", G = "G", T = "T", R = "[AG]", Y = "[CT]", S = "[GC]", W = "[AT]", K = "[GT]", M = "[AC]", B = "[CGT]", D = "[AGT]", H = "[ACT]", V = "[ACG]", N = "[ACGT]")
rx <- function(s) paste(iu[strsplit(s, "")[[1]]], collapse = "")
rc <- function(s) chartr("ACGT", "TGCA", vapply(strsplit(s, NULL), function(x) paste(rev(x), collapse = ""), ""))
readfa <- function(f) { l <- readLines(f); h <- grep("^>", l); setNames(l[h + 1], sub("^>", "", l[h])) }
sc <- fread(file.path(PD, "start_codons.tsv"))
core <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS"))
hits <- rbindlist(lapply(G, function(g) { s <- readfa(file.path(PD, sprintf("%s_1000.fa", g)))
  ok <- sc[genotype == g & start_codon == TRUE, gene]; s <- s[names(s) %in% ok]; r <- rc(s)
  rbindlist(lapply(seq_len(nrow(mot)), function(i) { p <- rx(mot$iupac[i])
    n <- lengths(regmatches(s, gregexpr(paste0("(?=", p, ")"), s, perl = TRUE))) + lengths(regmatches(r, gregexpr(paste0("(?=", p, ")"), r, perl = TRUE)))
    data.table(genotype = g, gene = names(s), motif = mot$motif[i], sites = n) })) }))
ft <- function(fg, d) { t <- table(factor(d$gene %in% fg, c(TRUE, FALSE)), factor(d$sites > 0, c(TRUE, FALSE))); f <- fisher.test(t)
  data.table(fg_genes = sum(d$gene %in% fg), fg_with = t[1, 1], bg_genes = sum(!d$gene %in% fg), bg_with = t[2, 1],
             odds_ratio = unname(f$estimate), ci_low = f$conf.int[1], ci_high = f$conf.int[2], p = f$p.value) }
# control (Col-0, Zhao CBF-activated, DRE/CRT)
es <- fread(ES); zact <- es[set == "zhao_cbf" & direction == "activated", unique(gene)]
d0 <- hits[genotype == "col0" & motif == "DRE_CRT" & !gene %in% core$col0]
ctl <- ft(zact, d0)[, `:=`(test = "Col-0 Zhao CBF-activated vs expressed background, DRE/CRT")]
z <- as.data.table(read_excel("data/external_gene_sets/PP2016-00533R1_Table_S3.xls", skip = 2, .name_repair = "minimal"))
z <- z[Type == "Activate", .(gene = toupper(`Gene ID`), zhao_sites = as.integer(`Number of CRT/DRE motif`))]
cmp <- merge(z, hits[genotype == "col0" & motif == "DRE_CRT", .(gene, our_sites = sites)], by = "gene")
ctl[, `:=`(zhao_genes_compared = nrow(cmp), zhao_share_with_site = round(mean(cmp$zhao_sites > 0), 3), our_share_with_site = round(mean(cmp$our_sites > 0), 3),
           presence_agreement = round(mean((cmp$zhao_sites > 0) == (cmp$our_sites > 0)), 3), count_agreement = round(mean(cmp$zhao_sites == cmp$our_sites), 3))]
fwrite(ctl, file.path(OUT, "cold_motif_control.tsv"), sep = "\t"); print(ctl)
if (!(ctl$p < 0.05 && ctl$odds_ratio > 1)) stop("positive control fails: the scan does not find DRE/CRT enriched in the CBF-activated genes")
# core vs background
res <- rbindlist(lapply(G, function(g) { fg <- core[[g]]
  rbindlist(lapply(mot$motif, function(m) ft(fg, hits[genotype == g & motif == m])[, `:=`(genotype = g, motif = m)])) }))
res[, q := p.adjust(p, "BH")][, call := fifelse(q < 0.05 & odds_ratio > 1, "enriched", fifelse(q < 0.05 & odds_ratio < 1, "depleted", "not enriched"))]
setcolorder(res, c("genotype", "motif")); fwrite(res, file.path(OUT, "cold_motif_enrichment.tsv"), sep = "\t")
print(res[, .(genotype, motif, fg = paste0(fg_with, "/", fg_genes), bg = paste0(bg_with, "/", bg_genes), OR = round(odds_ratio, 2),
              CI = sprintf("%.2f-%.2f", ci_low, ci_high), q = signif(q, 2), call)])
# CBF genes: symbol lookup in the TAIR alias file, orthogroup and core_status
al <- fread("data/annotation/gene_aliases_20140331.txt", header = TRUE, quote = "")
setnames(al, 1:2, c("locus", "symbol"))
cbf <- unique(al[toupper(symbol) %in% c("CBF1", "CBF2", "CBF3"), .(locus = toupper(locus), symbol)])
ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(get(ls(e)[1], e)) }
og <- unique(rbindlist(lapply(list.files("data/ComPlEx/PairwiseComparisons", "^comparison-col0-.*-table\\.RData$", full.names = TRUE),
  function(f) ld(f)[Species1 %in% cbf$locus, .(locus = Species1, OrthoGroup)])))
cs <- fread(CS)
cb <- merge(merge(cbf, og, by = "locus", all.x = TRUE), cs[, .(OrthoGroup, expressed_genes, all15, core_status)], by = "OrthoGroup", all.x = TRUE)
fwrite(cb, file.path(OUT, "cbf_genes_status.tsv"), sep = "\t", na = ""); print(cb)
