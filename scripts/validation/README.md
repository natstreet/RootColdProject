# Conserved-core validation

The tests and counts behind the conserved-core results. Each script is listed with its inputs and outputs in
`REPRODUCE.md`; run them from the repository root with the deposit unpacked into `data/`.

- `conservation_label_permutation.R` -- label-permutation null for the conserved core (Fig. S6a).
- `conservation_sizeaware_null.R` -- the same test with each orthogroup's family size kept (Results).
- `make_core_clusters.R`, then `string_external_validation.R` -- cluster coherence in STRING v12 (Fig. S6b);
  `string_by_cluster.R` breaks this down by cluster.
- `common_orthogroups.R`, `ecotype_deg_overlap.R`, `orthogroup_coexpressolog_comparisons.R` -- counts quoted in the
  Results.

Fig. S6 panel b needs the STRING v12 *A. thaliana* download (`3702.protein.links.detailed.v12.0.txt.gz` and
`3702.protein.aliases.v12.0.txt.gz`, https://string-db.org) in the working directory; it validates the clusters
against an independent database. With the four outputs present (`conservation_null_null.csv`,
`conservation_null_obs.txt`, `string_null.csv`, `string_obs.txt`), `Rscript scripts/figures/make_FigS6.R` draws
the figure.
