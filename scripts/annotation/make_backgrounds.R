#!/usr/bin/env Rscript
# make_backgrounds.R — per-species GO-enrichment backgrounds (the expressed-gene universe).
#
# A gene belongs to a species' background if it has any nonzero VST expression across the cold
# time-course samples (mean VST > 0); genes with all-zero expression are unexpressed and excluded.
# The same rule is applied to all six species.
#
# Input:  data/expression/<sp>_root_cold_stress_expression.txt   (VST expression, from the deposit)
# Output: data/annotation/bg_<sp>_<code>.rda                      (object bg_<code>, a character vector)
# Usage:  Rscript scripts/annotation/make_backgrounds.R

codes <- c(col0 = "ATC", ost0 = "ATO", aspen = "PT", birch = "BP", spruce = "PA", pine = "PS")
for (s in names(codes)) {
  ex <- read.delim(sprintf("data/expression/%s_root_cold_stress_expression.txt", s),
                   row.names = 1, check.names = FALSE)
  bg <- rownames(ex)[rowMeans(ex) > 0]
  obj <- paste0("bg_", codes[[s]]); assign(obj, bg)
  save(list = obj, file = sprintf("data/annotation/bg_%s_%s.rda", s, codes[[s]]))
  cat(sprintf("%-6s bg_%-3s : %d expressed genes\n", s, codes[[s]], length(bg)))
}
