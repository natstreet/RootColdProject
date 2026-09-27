#!/usr/bin/env Rscript
# check_figure_pngs.R — fails if any figure PNG has an alpha channel (RGBA or grey + alpha). Submission systems and PDF
# converters can flatten transparency to black, so every delivered PNG must be plain RGB (or grey) on white.
#
# Usage: Rscript scripts/figures/check_figure_pngs.R [file.png ...]   (default: every .png in . and figures/)

f <- commandArgs(trailingOnly = TRUE)
if (!length(f)) f <- c(Sys.glob("*.png"), Sys.glob("figures/*.png"))
if (!length(f)) stop("no PNG files to check")
bad <- character(0)
for (x in f) {
  d <- dim(png::readPNG(x))
  ch <- if (length(d) == 2) 1 else d[3]
  cat(sprintf("%-50s %4d x %4d  %s\n", basename(x), d[2], d[1], c("grey", "grey+alpha", "RGB", "RGBA")[ch]))
  if (ch %in% c(2, 4)) bad <- c(bad, x)
}
if (length(bad)) { cat("PNG with an alpha channel:", paste(basename(bad), collapse = ", "), "\n"); quit(status = 1) }
cat(sprintf("%d PNG files, none with an alpha channel\n", length(f)))
