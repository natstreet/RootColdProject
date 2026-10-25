# cluster_roles.R — manuscript cluster number of each core cluster's role (D095), from scripts/ComPlEx/cluster_numbering.tsv, checked
# against data/ComPlEx/cliques/manuscript_cluster_map.tsv (written from the same file by manuscript_cluster_map.R --numbering).
# Producers that test or report a particular cluster (a2 cell division vs translation, a1, cluster_profiles, manuscript_values) ask for
# it by role, so a renumbering changes their output names but never which cluster they test.
cluster_numbering <- function(file = "scripts/ComPlEx/cluster_numbering.tsv", map = "data/ComPlEx/cliques/manuscript_cluster_map.tsv") {
  nb <- utils::read.delim(file, stringsAsFactors = FALSE, colClasses = c(manuscript_cluster = "integer", analysis_cluster = "integer"))
  if (file.exists(map)) { mp <- utils::read.delim(map)
    if (!identical(mp$analysis_cluster[match(nb$manuscript_cluster, mp$manuscript_cluster)], nb$analysis_cluster))
      stop(map, " does not match ", file, "; run scripts/ComPlEx/manuscript_cluster_map.R --numbering") }
  nb
}
cluster_of <- function(role, nb = cluster_numbering()) {
  k <- nb$manuscript_cluster[nb$role == role]
  if (length(k) != 1) stop("no single cluster with role ", role, " in cluster_numbering.tsv")
  k
}
