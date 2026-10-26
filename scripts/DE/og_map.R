# og_map.R — gene -> orthogroup map per genotype for og_summary.R and make_Fig1.R (ledger R070, audit C003).
#   OG_MAP unset or "open" (D101): data/open/gene_orthogroup_map.tsv (the deposited map of Orthogroups_20240613, the table the Methods name).
#   OG_MAP = "gene2OG": data/ComPlEx/gene2OG/gene2OG_<genotype>.rda (carried over from the original analysis; it holds whole orthogroups
#     only, a subset whose rule is not documented; kept for comparison).
# Rows whose orthogroup is not an OrthoFinder ID (OG + digits) or whose gene is empty are dropped, and the result is checked, so no
# placeholder orthogroup ("No orthogroup", NA, "") is ever counted.
# og_map(s) returns data.frame(OrthoGroup, <s> = gene), one row per gene.
og_map <- function(s, which = Sys.getenv("OG_MAP", "open")) {
  if (which == "open") {
    sp <- c(col0 = "Arabidopsis_thaliana", ost0 = "Arabidopsis_thaliana", aspen = "Populus_tremula", birch = "Betula_pendula",
            spruce = "Picea_abies", pine = "Pinus_sylvestris")[[s]]
    if (!exists(".og_open", envir = globalenv())) assign(".og_open", utils::read.delim("data/open/gene_orthogroup_map.tsv", check.names = FALSE,
                                                         colClasses = "character")[, c("Ortholog_Group", "species", "gene")], envir = globalenv())
    m <- get(".og_open", envir = globalenv()); m <- unique(m[m$species == sp, c("Ortholog_Group", "gene")])
    m <- m[grepl("^OG[0-9]+$", m$Ortholog_Group) & !is.na(m$gene) & nzchar(m$gene), ]
    out <- stats::setNames(m, c("OrthoGroup", s))
  } else {
    fn <- c(col0 = "Col0", ost0 = "Ost0", aspen = "aspen", birch = "birch", spruce = "spruce", pine = "pine")[[s]]
    e <- new.env(); load(sprintf("data/ComPlEx/gene2OG/gene2OG_%s.rda", fn), envir = e); d <- as.data.frame(get(ls(e)[1], e))
    out <- stats::setNames(d[, 1:2], c("OrthoGroup", s))
  }
  stopifnot(all(grepl("^OG[0-9]+$", out$OrthoGroup)), !anyNA(out[[s]]), all(nzchar(out[[s]])))   # no placeholder orthogroup or gene
  out
}
