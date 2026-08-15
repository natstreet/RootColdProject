# Conserved-core validation and Figure S5

Two independent robustness tests of the conserved co-expression core, and the script that plots
them as Figure S5.

- `conservation_null_test.py` — degree-preserving permutation test (Figure S5a).
- `string_external_validation.py` — cluster coherence in STRING v12 (Figure S5b).
- `../figures/make_FigS5.py` — plots the two null distributions produced by the scripts above.

Unlike the other figures, Figure S5 is **not** regenerable from the figshare deposit with a single
command: panel (a) needs the per-species co-expression networks, and panel (b) needs an external
STRING download. Both are scripted; the steps are below.

## Panel (a) — `conservation_null_test.py`

Inputs it reads from the working directory:

- `networks/<species>.tsv` — one undirected edge list per species (`gene_a <tab> gene_b`), i.e. the
  mutual-rank co-expression networks. Build these with the ComPlEx pipeline in `scripts/ComPlEx/`
  (see the main README, "Reproducing the analysis", step 2) from the deposited `data/expression/`
  matrices, and write each species' edge list to `networks/<species>.tsv`.
- `orthogroups.tsv` — long-format orthogroup membership (`orthogroup <tab> species <tab> gene`),
  derived from the deposited wide table:

  ```r
  og <- read.delim("data/ComPlEx/Orthogroups_20240613.tsv", check.names = FALSE)
  sp <- c("Arabidopsis_thaliana","Betula_pendula","Populus_tremula","Picea_abies","Pinus_sylvestris")
  long <- do.call(rbind, lapply(sp, function(s)
    do.call(rbind, lapply(seq_len(nrow(og)), function(i) {
      g <- trimws(strsplit(og[i, s], ",")[[1]]); g <- g[nzchar(g)]
      if (length(g)) data.frame(orthogroup = og[i, 1], species = s, gene = g) })) ))
  write.table(long, "orthogroups.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  ```

Then: `python conservation_null_test.py` → writes `conservation_null_null.csv`,
`conservation_null_obs.txt`.

## Panel (b) — `string_external_validation.py`

Inputs it reads from the working directory:

- `3702.protein.links.detailed.v12.0.txt.gz` and `3702.protein.aliases.v12.0.txt.gz` — download once
  from STRING (https://string-db.org, organism 3702 = *A. thaliana*). External by design: the test
  validates the clusters against an independent database.
- `cliques.xlsx` — the conserved-core cliques with their Arabidopsis gene and cluster, exported from
  the deposit:

  ```r
  library(openxlsx)
  cc <- local({ e <- new.env(); load("data/ComPlEx/cliques/clique_HMcluster_list.rda", e); e[[ls(e)[1]]] })
  cl <- readRDS("data/ComPlEx/cliques/clique_genes_filterable_COMPLETE.RDS")
  col0 <- tapply(c(cl$GeneSpecies1, cl$GeneSpecies2), c(cl$OrthoGroup, cl$OrthoGroup),
                 function(g) g[grepl("^AT", g)][1])   # one Col-0 gene per orthogroup
  cc$OG <- sub("_.*$", "", cc$ID); cc$Col0 <- col0[cc$OG]
  write.xlsx(cc[, c("Col0","Cluster")], "cliques.xlsx")
  ```

Then: `python string_external_validation.py` → writes `string_null.csv`, `string_obs.txt`.

## Figure S5

With the four output files present:

```
python ../figures/make_FigS5.py     # writes Figure_S5.pdf / .png
```
