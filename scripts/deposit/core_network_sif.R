#!/usr/bin/env Rscript
# core_network_sif.R — conserved_core_network.sif: the conserved core as a cross-species graph for Cytoscape. Edges are
# the co-expressolog pairs (co_expressologs.tsv) of the core orthogroups (conserved_core.tsv); nodes are species:gene.
# Columns: source, target, OrthologGroup, MaxpVal, effectSize, in the row order of co_expressologs.tsv.
# Usage: Rscript scripts/deposit/core_network_sif.R <deposit data/open dir>   (reads and writes there)
suppressMessages(library(data.table))
D <- commandArgs(trailingOnly = TRUE)[1]
ce <- fread(file.path(D, "co_expressologs.tsv"), colClasses = c(MaxpVal = "character", effectSize = "character"))
core <- fread(file.path(D, "conserved_core.tsv"))
x <- ce[OrthologGroup %in% core$OrthoGroup]
sif <- x[, .(source = paste0(Species1, ":", GeneSpecies1), target = paste0(Species2, ":", GeneSpecies2), OrthologGroup, MaxpVal, effectSize)]
if (!setequal(unique(sif$OrthologGroup), unique(core$OrthoGroup))) stop("a core orthogroup has no co-expressolog pair")
fwrite(sif, file.path(D, "conserved_core_network.sif"), sep = "\t", quote = FALSE)
cat(sprintf("conserved_core_network.sif: %d edges over %d nodes, %d orthogroups\n", nrow(sif),
            uniqueN(c(sif$source, sif$target)), uniqueN(sif$OrthologGroup)))
