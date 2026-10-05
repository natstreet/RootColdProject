#!/usr/bin/env Rscript
# orthogroup_size_flag.R — how likely an orthogroup of a given size is to be conserved in all 15 comparisons, and that
# probability for each core orthogroup (NRS, 2026-10-05: the size flag replaces the excluded large orthogroups).
#
# Expressed genes per orthogroup: the genes of the orthogroup table (Arabidopsis, birch, Norway spruce, Scots pine,
# aspen) that are rows of the six expression matrices; an Arabidopsis gene counts once whichever ecotype it is in.
# Eligible orthogroups: at least one expressed gene in each of the five species.
# Response: a co-expressolog (Max.p.val < 0.1) in every one of the 15 genotype comparisons.
# Model: logistic regression of the response on log(expressed genes), over the eligible orthogroups.
#
# Inputs:  data/ComPlEx/Orthogroups_20240613.tsv, data/expression/<genotype>_root_cold_stress_expression.txt,
#          data/ComPlEx/all_expressed_genes.RData, data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS,
#          data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS
# Outputs: results/orthogroup_size_model.tsv   coefficients, n, and the fitted probability at fixed sizes
#          results/orthogroup_size_flag.tsv    per core orthogroup: expressed_genes, complete_cliques, p_all15_given_size
#          results/orthogroup_sizes_all.tsv    per orthogroup: expressed genes per species, eligible, all15
# Usage (repository root): Rscript scripts/validation/orthogroup_size_flag.R
suppressPackageStartupMessages(library(dplyr))
dir.create("results", showWarnings = FALSE)
P_COEXP <- 0.1
SIZES <- c(10, 20, 50, 100, 200, 300)
GENOTYPES <- c("col0", "ost0", "aspen", "birch", "spruce", "pine")
SPECIES <- c(Arabidopsis_thaliana = "arabidopsis", Betula_pendula = "birch", Picea_abies = "spruce",
             Pinus_sylvestris = "pine", Populus_tremula = "aspen")

expressed <- unique(unlist(lapply(GENOTYPES, function(g)
  read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", g))[[1]])))
og <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE, colClasses = "character")
names(og)[1] <- "OrthoGroup"
if (anyDuplicated(og$OrthoGroup)) stop("duplicated orthogroups")
long <- do.call(rbind, lapply(names(SPECIES), function(s) {
  g <- strsplit(og[[s]], ", ", fixed = TRUE)
  data.frame(OrthoGroup = rep(og$OrthoGroup, lengths(g)), species = SPECIES[[s]], gene = unlist(g))
})) %>% filter(gene != "", gene %in% expressed) %>% distinct()
per_sp <- long %>% count(OrthoGroup, species) %>% tidyr::pivot_wider(names_from = species, values_from = n, values_fill = 0)
for (s in SPECIES) if (!s %in% names(per_sp)) per_sp[[s]] <- 0L
sizes <- per_sp %>% mutate(expressed_genes = rowSums(across(all_of(unname(SPECIES)))),
                           eligible = if_all(all_of(unname(SPECIES)), ~ .x > 0))

load("data/ComPlEx/all_expressed_genes.RData")                                       # expr_genes
pairs <- expr_genes %>% filter(MaxpVal < P_COEXP) %>%
  mutate(cmp = ifelse(Species1 < Species2, paste(Species1, Species2), paste(Species2, Species1))) %>%
  distinct(OrthologGroup, cmp) %>% count(OrthologGroup, name = "comparisons")
n_cmp <- choose(length(GENOTYPES), 2)
sizes <- sizes %>% left_join(pairs, by = c(OrthoGroup = "OrthologGroup")) %>%
  mutate(comparisons = coalesce(comparisons, 0L), all15 = comparisons == n_cmp)

d <- filter(sizes, eligible)
fit <- glm(all15 ~ log(expressed_genes), family = binomial, data = d)
if (!fit$converged) stop("logistic regression did not converge")
b <- coef(fit)
model <- rbind(
  data.frame(quantity = c("n_eligible", "n_all15", "intercept", "slope_log_expressed_genes"),
             value = c(nrow(d), sum(d$all15), b[[1]], b[[2]])),
  data.frame(quantity = paste0("p_all15_at_", SIZES),
             value = unname(predict(fit, data.frame(expressed_genes = SIZES), type = "response"))))
write.table(model, "results/orthogroup_size_model.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sizes, "results/orthogroup_sizes_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

core <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_filtered.RDS")$OrthoGroup
cliques <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS") %>% count(OrthoGroup, name = "complete_cliques")
flag <- data.frame(OrthoGroup = sort(unique(core))) %>%
  left_join(select(sizes, OrthoGroup, expressed_genes, eligible), by = "OrthoGroup") %>%
  left_join(cliques, by = "OrthoGroup")
if (anyNA(flag$expressed_genes) || !all(flag$eligible)) stop("a core orthogroup lacks an expressed gene in some species")
flag$p_all15_given_size <- unname(predict(fit, flag, type = "response"))
write.table(select(flag, -eligible), "results/orthogroup_size_flag.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("eligible %d orthogroups (%d conserved in all 15); intercept %.4f, slope %.4f; core %d orthogroups, expressed genes %d-%d\n",
            nrow(d), sum(d$all15), b[[1]], b[[2]], nrow(flag), min(flag$expressed_genes), max(flag$expressed_genes)))
print(model, row.names = FALSE)
