#!/usr/bin/env Rscript
# external_gene_sets.R — Arabidopsis cold gene sets from published supplementary tables (DECISIONS.md D062 item 3, D063;
# design R040 as amended in A035, N131). Files as downloaded (ledger/attachments/external_gene_sets/MANIFEST.tsv); each
# set's size is checked against the number its paper states, and the script stops on a mismatch.
#   jia_cbf       Jia et al. (2016) New Phytol. 10.1111/nph.14088, Table S7: CBF-regulated genes, 112 activated + 22
#                 repressed = 134 (abstract)
#   zhao_cbf      Zhao et al. (2016) Plant Physiol. 10.1104/pp.16.00533, Table S3: 414 CBF-regulated cold-responsive genes,
#                 346 activated + 68 repressed (abstract)
#   cbf_both      jia_cbf and zhao_cbf
#   calixto_up / calixto_down  Calixto et al. (2018) Plant Cell 10.1105/tpc.18.00177, Supplemental Data Set 2: genes up
#                 (down) in at least one contrast group; Data Set 1B All_DE_genes = 7302 (text) is checked as their union
# Menges et al. (2003) is not built: its 1,082-gene set (Table 2) is not among the supplementary files and cannot be
# reproduced from Tables S5/S11 (ledger N135, N145).
# Output (<out dir>): external_gene_sets.tsv (set, gene, direction), external_gene_sets_counts.tsv
# Usage: Rscript external_gene_sets.R <attachments dir> <dir with Zhao Table S3 .xls and Calixto Data Sets 1, 2 .xlsx> <out dir>
suppressPackageStartupMessages({ library(data.table); library(readxl) })
a <- commandArgs(trailingOnly = TRUE); ATT <- a[1]; RAW <- a[2]; OUT <- a[3]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
isat <- function(v) grepl("^AT[1-5MC]G[0-9]{5}$", toupper(v))
chk <- function(lab, got, want) { cat(sprintf("%-28s %6d (paper %d)\n", lab, got, want)); if (got != want) stop(lab, ": count differs from the paper") }
j <- read_excel(file.path(ATT, "nph14088-sup-0007-tables7.xlsx"), col_names = FALSE, .name_repair = "minimal")[[1]]
rep_row <- which(j == "Repressed")
jia <- data.table(gene = toupper(j), direction = ifelse(seq_along(j) > rep_row, "repressed", "activated"))[isat(gene)]
chk("Jia activated", jia[direction == "activated", uniqueN(gene)], 112); chk("Jia repressed", jia[direction == "repressed", uniqueN(gene)], 22)
z <- as.data.table(read_excel(file.path(RAW, "PP2016-00533R1_Table_S3.xls"), skip = 2, .name_repair = "minimal"))
zhao <- z[isat(`Gene ID`), .(gene = toupper(`Gene ID`), direction = ifelse(Type == "Activate", "activated", ifelse(Type == "Repress", "repressed", NA)))]
chk("Zhao activated", zhao[direction == "activated", uniqueN(gene)], 346); chk("Zhao repressed", zhao[direction == "repressed", uniqueN(gene)], 68)
d2 <- read_excel(file.path(RAW, "Supplemental_Dataset_2.xlsx"), col_names = FALSE, .name_repair = "minimal")
hdr <- as.character(unlist(d2[6, ])); body <- d2[-(1:6), ]
cal <- rbindlist(lapply(seq_along(hdr), function(i) data.table(gene = toupper(body[[i]]), direction = hdr[i])))[isat(gene) & direction %in% c("up", "down")]
d1 <- read_excel(file.path(RAW, "Supplemental_Dataset_1.xlsx"), sheet = "Supplemental Dataset 1B", skip = 4)
alld <- unique(toupper(na.omit(d1$All_DE_genes)))
chk("Calixto All_DE_genes", length(alld), 7302)
up <- unique(cal[direction == "up", gene]); dn <- unique(cal[direction == "down", gene])
cat(sprintf("Calixto up %d, down %d, both %d; union %d, of them in All_DE_genes %d\n", length(up), length(dn), length(intersect(up, dn)),
            length(union(up, dn)), length(intersect(union(up, dn), alld))))
sets <- rbind(jia[, .(set = "jia_cbf", gene, direction)], zhao[, .(set = "zhao_cbf", gene, direction)],
              data.table(set = "cbf_both", gene = intersect(jia$gene, zhao$gene), direction = NA_character_),
              data.table(set = "calixto_up", gene = up, direction = "up"), data.table(set = "calixto_down", gene = dn, direction = "down"),
              data.table(set = "calixto_all_de", gene = alld, direction = NA_character_))
fwrite(sets, file.path(OUT, "external_gene_sets.tsv"), sep = "\t")
cnt <- sets[, .(genes = uniqueN(gene)), by = set]; fwrite(cnt, file.path(OUT, "external_gene_sets_counts.tsv"), sep = "\t"); print(cnt, row.names = FALSE)
