#!/usr/bin/env Rscript
# core_module_directions.R — the conserved co-expression neighbours of each conserved-core clique member, and whether
# a member whose cold response differs in direction between genotypes changes direction together with them.
#
# Conserved neighbours: ComPlEx tests each orthologue pair by the overlap of the two genes' network neighbourhoods
# (scripts/ComPlEx/core_neighbourhoods.R recovers the overlapping genes). A neighbour of a member in one genotype is
# conserved when it lies in that overlap in at least MIN_PAIRS of the member's five pairwise comparisons.
# ComPlEx uses |r|, so a conserved neighbour can be positively or negatively correlated with the member. For each
# member, the Pearson correlation with each conserved neighbour is taken over all samples of its genotype
# (data/expression/<sp>_root_cold_stress_expression.txt, the ComPlEx input).
#
# Outputs (tables/):
#   core_conserved_neighbours.tsv.gz  OrthoGroup, species, member, neighbour, n_pairs, min_pairs (the MIN_PAIRS used),
#                                     r (member-neighbour correlation)
#   core_module_directions.tsv        one row per clique and genotype: number of conserved neighbours, the fraction
#                                     positively correlated with the member, the member's response (from
#                                     tables/core_response_patterns.tsv), and the median log2FC of the positively and
#                                     of the negatively correlated neighbours at each time point
#   core_module_sign_agreement.tsv    one row per clique: whether a neighbour orthogroup keeps the same sign of
#                                     correlation with the member across genotypes. For every orthogroup with conserved
#                                     neighbours in two or more genotypes, its sign in a genotype is that of the mean r of
#                                     its conserved neighbours there; agreement is the fraction of genotype pairs (over
#                                     all such orthogroups) in which the two signs agree. High agreement in a clique whose
#                                     members respond in opposite directions means the neighbours change direction with
#                                     the member; low agreement means the member's relationship to them has changed.
#                                     chance_agreement is the agreement expected if signs were independent between
#                                     genotypes, given each genotype's share of positive orthogroups; excess = observed
#                                     minus chance.
#                                     Orthogroups from data/ComPlEx/orthologs-<sp1>-<sp2>.RData (the ComPlEx orthology).
# Usage:  Rscript scripts/analysis/core_module_directions.R [MIN_PAIRS, default 3]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); MIN_PAIRS <- if (length(a)) as.integer(a[1]) else 3L
SP <- c("col0", "ost0", "aspen", "birch", "spruce", "pine"); TPS <- c("6h", "24h", "3d", "10d")
ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

nb <- rbindlist(lapply(list.files("data/ComPlEx/core_neighbourhoods", "[.]tsv[.]gz$", full.names = TRUE), fread))
if (uniqueN(nb[, .(OrthoGroup, species)]) != 132 * 6 || nb[, uniqueN(species), by = .(OrthoGroup)][, any(V1 != 6)])
  stop("core_neighbourhoods: expected all 15 species pairs for every core clique")
cn <- nb[overlap == TRUE, .(n_pairs = .N), by = .(OrthoGroup, species, member, neighbour)][n_pairs >= MIN_PAIRS]
cn[, min_pairs := MIN_PAIRS]

expr <- lapply(setNames(SP, SP), function(s) {
  x <- fread(sprintf("data/expression/%s_root_cold_stress_expression.txt", s))
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m })
cn[, r := { e <- expr[[species[1]]]; as.vector(cor(t(e[neighbour, , drop = FALSE]), e[member[1], ])) },
   by = .(OrthoGroup, species, member)]
dir.create("tables", showWarnings = FALSE)
fwrite(cn, "tables/core_conserved_neighbours.tsv.gz", sep = "\t")

de <- fread("DE/de_results.tsv")[timepoint %in% TPS]
pat <- fread("tables/core_response_patterns.tsv")
lfc <- dcast(de, species + gene ~ timepoint, value.var = "log2FoldChange")
cl <- merge(cn, lfc, by.x = c("species", "neighbour"), by.y = c("species", "gene"), all.x = TRUE)
med <- function(v) if (length(v)) round(median(v, na.rm = TRUE), 3) else NA_real_
out <- cl[, c(list(n_conserved = .N, frac_positive = round(mean(r > 0), 3)),
              setNames(lapply(TPS, function(t) med(get(t)[r > 0])), paste0("pos_lfc_", TPS)),
              setNames(lapply(TPS, function(t) med(get(t)[r < 0])), paste0("neg_lfc_", TPS))),
          by = .(OrthoGroup, species, member)]
out <- merge(CJ(OrthoGroup = unique(nb$OrthoGroup), species = SP), out, by = c("OrthoGroup", "species"), all.x = TRUE)
out[is.na(n_conserved), n_conserved := 0L]
long <- melt(pat[, c("OrthoGroup", paste0("response_", SP)), with = FALSE], id.vars = "OrthoGroup",
             variable.name = "species", value.name = "member_response")
long[, species := sub("^response_", "", species)]
out <- merge(out, long, by = c("OrthoGroup", "species"))
out[, species := factor(species, SP)]; setorder(out, OrthoGroup, species)
fwrite(out, "tables/core_module_directions.tsv", sep = "\t")
## ---- sign agreement of neighbour orthogroups across genotypes ----
og <- unique(rbindlist(lapply(list.files("data/ComPlEx", "^orthologs-.*[.]RData$", full.names = TRUE), function(f) {
  o <- ld(f); rbind(data.table(gene = o$Species1, og = o$OrthoGroup), data.table(gene = o$Species2, og = o$OrthoGroup)) })))
if (anyDuplicated(og$gene)) stop("a gene in two orthogroups")
cn[, og := og$og[match(neighbour, og$gene)]]
if (anyNA(cn$og)) stop("conserved neighbours without an orthogroup")
s_og <- cn[, .(sgn = sign(mean(r))), by = .(OrthoGroup, og, species)]
s_og[, p_pos := mean(sgn > 0), by = .(OrthoGroup, species)]     # the genotype's share of positive orthogroups
agree <- s_og[, { if (.N < 2) NULL else { p <- combn(.N, 2)
  list(same = sgn[p[1, ]] == sgn[p[2, ]],
       chance = p_pos[p[1, ]] * p_pos[p[2, ]] + (1 - p_pos[p[1, ]]) * (1 - p_pos[p[2, ]])) } }, by = .(OrthoGroup, og)]
ag <- agree[, .(neighbour_orthogroups = uniqueN(og), genotype_pairs = .N, sign_agreement = round(mean(same), 3),
                chance_agreement = round(mean(chance), 3)), by = OrthoGroup]
ag[, excess := sign_agreement - chance_agreement]
ag <- merge(ag, pat[, .(OrthoGroup, n_genotypes_DE, lineage, direction)], by = "OrthoGroup")
setorder(ag, -n_genotypes_DE, OrthoGroup)
fwrite(ag, "tables/core_module_sign_agreement.tsv", sep = "\t")
cat("Sign agreement of neighbour orthogroups across genotypes, by the members' response direction (2+ genotypes DE):\n")
print(ag[n_genotypes_DE >= 2, .(cliques = .N, median_agreement = median(sign_agreement),
                                median_chance = median(chance_agreement), median_excess = median(excess),
                                cliques_below_chance = sum(excess <= 0)), by = direction])
print(ag[OrthoGroup == "OG0000166"])
cat(sprintf("all cliques: median agreement %.3f (range %.3f-%.3f)\n", median(ag$sign_agreement), min(ag$sign_agreement), max(ag$sign_agreement)))

cat(sprintf("conserved neighbours (overlap in >= %d of 5 comparisons): median %s per member (range %d-%d); %d members with none\n",
            MIN_PAIRS, median(out$n_conserved), min(out$n_conserved), max(out$n_conserved), sum(out$n_conserved == 0)))
cat("wrote tables/core_conserved_neighbours.tsv.gz, tables/core_module_directions.tsv, tables/core_module_sign_agreement.tsv\n")
