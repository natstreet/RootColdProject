#!/usr/bin/env Rscript
# menges2005_core_regulators.R — descriptive (no test; ledger N131): which of the core cell-cycle regulators of Menges
# et al. (2005) Plant J. 10.1111/j.1365-313x.2004.02319.x (Table S1, sheet "CC genes_all expression": MIPS ID, synonym)
# are Col-0 members of the conserved core's complete cliques, with orthogroup, cluster, stability status and whether the
# gene is in the representative clique.
# Output (<out dir>): menges2005_core_regulators_in_core.tsv. Usage (working tree): Rscript menges2005_core_regulators.R <tpj_2319_sm_tables1.xls> <kept_lost_orthogroups.tsv> <out dir>
suppressPackageStartupMessages({ library(data.table); library(readxl) })
a <- commandArgs(trailingOnly = TRUE); dir.create(a[3], recursive = TRUE, showWarnings = FALSE)
s1 <- as.data.table(read_excel(a[1], sheet = "CC genes_all expression", .name_repair = "minimal"))
reg <- unique(s1[grepl("^AT[1-5MC]G[0-9]{5}$", toupper(`MIPS ID`)), .(gene = toupper(`MIPS ID`), synonym)])
cat("Menges 2005 Table S1 regulators:", nrow(reg), "\n")
kl <- fread(a[2]); w <- as.data.table(readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS"))[OrthoGroup %in% kl$OrthoGroup]
m <- unique(w[col0 %in% reg$gene, .(gene = col0, OrthoGroup)])
m <- merge(m, reg, by = "gene")
m <- merge(m, kl[, .(OrthoGroup, cluster, core_status, representative_col0 = col0_gene)], by = "OrthoGroup")
m[, in_representative := gene == representative_col0][, representative_col0 := NULL]
setorder(m, cluster, core_status, synonym); fwrite(m, file.path(a[3], "menges2005_core_regulators_in_core.tsv"), sep = "\t"); print(m, row.names = FALSE)
