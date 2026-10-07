# loo_task.R — one ComPlEx run of the D048/D052 leave-one-sample-out test (ledger N115, R039, D059). Run inside an empty task
# directory: links the inputs, writes the left-out genotype's matrix without the sample (none for a baseline), runs the
# analysis's ComPlEx.R unchanged, and keeps the comparison rows of the orthogroups with a complete clique.
# Usage: Rscript task.R <work dir> <task id>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); W <- a[1]; id <- as.integer(a[2])
tk <- fread(file.path(W, "tasks.tsv"))[task == id]; I <- file.path(W, "inputs/data")
dir.create("data/ComPlEx", recursive = TRUE); dir.create("data/expression")
for (p in c("orthologs", "annotation")) file.symlink(file.path(I, p), file.path("data", p))
of <- sprintf("orthologs-%s-%s.RData", tk$s1, tk$s2); file.symlink(file.path(I, "ComPlEx", of), file.path("data/ComPlEx", of))
for (s in c(tk$s1, tk$s2)) {
  f <- sprintf("%s_root_cold_stress_expression.txt", s)
  if (s == tk$left_out_genotype) {
    e <- fread(file.path(I, "expression", f)); if (!tk$sample %in% names(e)) stop("sample not in matrix")
    e[, (tk$sample) := NULL]; fwrite(e, file.path("data/expression", f), sep = "\t", quote = FALSE)
  } else file.symlink(file.path(I, "expression", f), file.path("data/expression", f))
}
st <- system(sprintf("Rscript %s %s %s > complex.log 2>&1", file.path(W, "ComPlEx.R"), tk$s1, tk$s2))
tf <- sprintf("data/ComPlEx/PairwiseComparisons/comparison-%s-%s-abspearsonMR0.03-table.RData", tk$s1, tk$s2)
if (st != 0 || !file.exists(tf)) { system("tail -20 complex.log"); stop("ComPlEx failed") }
e <- new.env(); load(tf, envir = e); x <- as.data.table(get(ls(e)[1], e))
ogs <- readLines(file.path(W, "clique_ogs.txt"))
x <- x[OrthoGroup %in% ogs, .(Species1, Species2, OrthoGroup, Species1.neigh.overlap, Species2.neigh.overlap, Max.p.val)]
fwrite(x, file.path(W, "out", sprintf("task_%03d.tsv.gz", id)), sep = "\t")
cat(sprintf("task %d: %s-%s, left out %s %s: %d rows kept\n", id, tk$s1, tk$s2, tk$left_out_genotype, tk$sample, nrow(x)))
