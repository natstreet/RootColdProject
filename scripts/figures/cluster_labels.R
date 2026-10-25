# cluster_labels.R — approved display names of core clusters (scripts/figures/cluster_labels_approved.tsv; one row per named cluster,
# with the decision that named it). Used by make_Fig4.R (via fig4_captions_family_level.R), make_FigS5.R, make_Fig5.R and the Table S3
# legend. approved_labels() returns the table; label_for() gives the approved name of a cluster when the figure names the same GO term,
# the Table S3 term name otherwise, and stops when a figure names a different GO term for an approved cluster.
approved_labels <- function(file = "scripts/figures/cluster_labels_approved.tsv")
  if (file.exists(file)) utils::read.delim(file, stringsAsFactors = FALSE, colClasses = c(cluster = "integer")) else
    data.frame(cluster = integer(0), go_id = character(0), label = character(0), decision = character(0))
label_for <- function(cluster, go_id, term, ap = approved_labels()) {
  i <- match(cluster, ap$cluster)
  if (is.na(i)) return(paste0(toupper(substr(term, 1, 1)), substring(term, 2)))
  if (ap$go_id[i] != go_id) stop("cluster ", cluster, ": figure names ", go_id, ", approved label is for ", ap$go_id[i], " (", ap$decision[i], ")")
  ap$label[i]
}
