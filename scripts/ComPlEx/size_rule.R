#!/usr/bin/env Rscript
# size_rule.R — the even-odds orthogroup size rule (DECISIONS.md D026, D043): an orthogroup with a complete clique is
# not counted as conserved when its size alone makes a complete clique more likely than not. Run after
# complete_cliques.R and before core_cliques.R, which takes this table as its second argument.
#
# Expressed genes per orthogroup: the genes of the orthogroup table (Arabidopsis, birch, Norway spruce, Scots pine,
# aspen) that are rows of the six expression matrices; an Arabidopsis gene counts once whichever ecotype it is in (the
# definition of scripts/validation/orthogroup_size_flag.R).
# Eligible orthogroups: at least one expressed gene in each of the five species.
# Response: the orthogroup has a complete clique (complete_cliques.R: one gene per genotype, all 15 pairs co-expressologs).
# Model: logistic regression of the response on log(expressed genes), over the eligible orthogroups. Even odds is the
# size where the fitted probability is 0.5, exp(-intercept / slope); its 95% interval is the percentile interval over
# bootstrap resamples of the eligible orthogroups.
# Rule: counted_as_conserved = has a complete clique and fitted probability < 0.5. Optional second argument: a fixed size
# threshold in expressed genes instead (sensitivity runs); an orthogroup at or above it is not counted.
#
# Inputs:  data/ComPlEx/Orthogroups_20240613.tsv, data/expression/<genotype>_root_cold_stress_expression.txt,
#          data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS
# Outputs: <out>.tsv        one row per eligible orthogroup: OrthoGroup, expressed_genes, has_complete_clique,
#                           complete_cliques, p_core_given_size, counted_as_conserved
#          <out>_model.tsv  n, coefficients, even-odds size, bootstrap interval, resamples, seed, the orthogroups not counted
# Usage (repository root): Rscript scripts/ComPlEx/size_rule.R [out prefix, default results/orthogroup_size_rule] [threshold]
suppressPackageStartupMessages(library(dplyr))
a <- commandArgs(trailingOnly = TRUE)
out <- if (length(a) >= 1) a[1] else "results/orthogroup_size_rule"
threshold <- if (length(a) >= 2) as.numeric(a[2]) else NA
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
N_BOOT <- 2000; SEED <- 1
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

cl <- readRDS("data/ComPlEx/cliques/conserved_genes_COMPLETE_wider.RDS") %>% count(OrthoGroup, name = "complete_cliques")
d <- sizes %>% filter(eligible) %>% select(OrthoGroup, expressed_genes) %>%
  left_join(cl, by = "OrthoGroup") %>%
  mutate(complete_cliques = coalesce(complete_cliques, 0L), has_complete_clique = complete_cliques > 0)
if (!all(cl$OrthoGroup %in% d$OrthoGroup)) stop("an orthogroup with a complete clique is not eligible")

fit <- glm(has_complete_clique ~ log(expressed_genes), family = binomial, data = d)
if (!fit$converged) stop("logistic regression did not converge")
b <- coef(fit); even <- exp(-b[[1]] / b[[2]])
set.seed(SEED)
boot <- replicate(N_BOOT, {
  s <- d[sample.int(nrow(d), replace = TRUE), ]
  f <- suppressWarnings(glm(has_complete_clique ~ log(expressed_genes), family = binomial, data = s))
  exp(-coef(f)[[1]] / coef(f)[[2]])
})
ci <- unname(quantile(boot, c(0.025, 0.975)))

d$p_core_given_size <- unname(predict(fit, d, type = "response"))
d$counted_as_conserved <- d$has_complete_clique &
  (if (is.na(threshold)) d$p_core_given_size < 0.5 else d$expressed_genes < threshold)
d <- d %>% select(OrthoGroup, expressed_genes, has_complete_clique, complete_cliques, p_core_given_size, counted_as_conserved) %>%
  arrange(desc(expressed_genes), OrthoGroup)
write.table(d, paste0(out, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

flagged <- d %>% filter(has_complete_clique, !counted_as_conserved)
model <- data.frame(
  quantity = c("n_eligible", "n_with_complete_clique", "intercept", "slope_log_expressed_genes", "even_odds_expressed_genes",
               "bootstrap_ci_low", "bootstrap_ci_high", "bootstrap_resamples", "seed", "rule",
               "n_counted_as_conserved", "not_counted"),
  value = c(nrow(d), sum(d$has_complete_clique), b[[1]], b[[2]], even, ci[1], ci[2], N_BOOT, SEED,
            if (is.na(threshold)) "fitted probability < 0.5" else paste("expressed genes <", threshold),
            sum(d$counted_as_conserved), paste(flagged$OrthoGroup, collapse = ",")))
write.table(model, paste0(out, "_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(model, row.names = FALSE)
print(as.data.frame(flagged), row.names = FALSE)
