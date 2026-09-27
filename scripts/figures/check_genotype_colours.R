#!/usr/bin/env Rscript
# check_genotype_colours.R — fails if any script other than genotype_colours.R assigns a colour code to a genotype.
#
# Looks for a genotype name or label (col0, Col-0, ost0, Ost-0, aspen, birch, spruce, pine, in any case) assigned a
# hex colour, e.g. `col0 = "#0072B2"` or `"Pine"="#e377c2"`, in every R, Python and HTML file under scripts/, and
# checks that scripts/figures/genotype_colours.json matches genotype_colours.R.
#
# Usage: Rscript scripts/figures/check_genotype_colours.R      (exit status 1 on any finding)

files <- list.files("scripts", pattern = "[.](R|Rmd|py|html)$", recursive = TRUE, full.names = TRUE)
files <- files[basename(files) != "genotype_colours.R" & basename(files) != "check_genotype_colours.R"]
pat <- "[\"'`]?(col-?0|ost-?0|aspen|birch|spruce|pine)[\"'`]?\\s*[=:]\\s*[\"']#[0-9A-Fa-f]{6}"
bad <- character(0)
for (f in files) {
  l <- readLines(f, warn = FALSE); h <- grep(pat, l, ignore.case = TRUE, perl = TRUE)
  if (length(h)) bad <- c(bad, sprintf("%s:%d: %s", f, h, trimws(l[h])))
}
source("scripts/figures/genotype_colours.R")
js <- readLines("scripts/figures/genotype_colours.json", warn = FALSE)
got <- regmatches(js, regexec('"([a-z0-9]+)": "(#[0-9A-Fa-f]{6})"', js))
got <- do.call(rbind, got[lengths(got) == 3])
json_ok <- !is.null(got) && identical(setNames(got[, 3], got[, 2]), GENOTYPE_COLOURS)
if (length(bad)) cat("genotype colours defined outside genotype_colours.R:\n", paste(" ", bad, collapse = "\n"), "\n")
if (!json_ok) cat("genotype_colours.json does not match genotype_colours.R; rerun Rscript scripts/figures/genotype_colours.R\n")
if (length(bad) || !json_ok) quit(status = 1)
cat(sprintf("genotype colours: %d files checked, none defines its own; JSON matches\n", length(files)))
