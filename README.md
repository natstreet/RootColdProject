# RootColdProject

Analysis code for:

**A comparative transcriptomic resource for the root cold response of boreal trees across the angiosperm–gymnosperm divide**
Norway spruce (*Picea abies*), Scots pine (*Pinus sylvestris*), silver birch (*Betula pendula*),
aspen (*Populus tremula*) and two *Arabidopsis thaliana* ecotypes (Col-0, Ost-0).

This repository contains the scripts needed to reproduce the analyses and figures in the paper.
The underlying data are archived separately (see **Data**); download them into a local `data/`
directory before running the scripts.

## Overview
Fine-root transcriptomes of four boreal trees and two *Arabidopsis* ecotypes were profiled over a
ten-day 5 °C (chilling) treatment. Differential expression is largely lineage-specific, but
comparative co-expression (ComPlEx) recovers a conserved regulatory core of 256 co-expressolog
orthogroups (116 of which form complete cliques) shared across all 15 pairwise species
comparisons and spanning more than 300 million years of evolution. Supporting analyses test that this core is significant against a
permutation null, is coherent in an independent co-expression compendium (STRING), replicates in
independent data, and shows conserved promoter cis-motifs.

## Data
The scripts read their inputs from a local `data/` directory. Most of it is not held in this
repository and is obtained as follows:

- **Processed data** (per-species VST expression matrices, sample metadata, DEG tables,
  transcription-factor lists, ComPlEx outputs and clique definitions, functional annotation,
  the OrthoFinder orthogroup table): SciLifeLab figshare, https://figshare.scilifelab.se/
  (DOI:10.17044/scilifelab.32747586). Download `data.tar.gz` and unpack it; it expands to a
  `data/` directory containing `expression/`, `metadata/`, `DEGs/`, `TranscriptionFactors/`,
  `annotation/`, `orthologs/` and `ComPlEx/`.
  The `ComPlEx/` step can be re-run either from the cached per-pair orthologue
  tables (`ComPlEx/orthologs-*.RData`) or, from scratch, from `orthologs/orthologs.parquet`
  together with `annotation/gene_aliases_20140331.txt`. *Arabidopsis* annotation is derived
  from `orthologs.parquet` itself, so that it cannot disagree with the orthology used for the
  comparisons. The variance-stabilised
  `expression/` matrices are the starting point for network construction; the per-species
  `DESeqDataSet` objects (`dds/`) are included and are the input to
  `scripts/DE/differential_expression.R`.
- **Raw sequencing reads:** ENA umbrella BioProject **PRJEB104158** (components PRJEB104117
  Arabidopsis, PRJEB104118 birch, PRJEB104119 aspen, PRJEB104120 Scots pine, PRJEB26918 Norway
  spruce roots).
- **Genome assemblies and annotation releases:** *Picea abies* v2.0 and *Pinus sylvestris* v1.0
  (Kalman et al. 2025), *Populus tremula* v2.2, *Betula pendula* v1.4, and for *A. thaliana* the
  Araport11 June 2016 release (Cheng et al. 2017), quantified against the `Araport11_all.201606.cdna`
  transcript fasta. The Arabidopsis release was established by matching the transcript identifiers in
  the deposited Salmon output against candidate references: the 53,467 quantified transcripts are
  exactly that fasta's 54,013 minus 546 exact-sequence duplicates collapsed at index time. Gene
  identifiers in the deposited tables belong to these releases; joining them to a resource built on a
  different release will silently fail to match rather than error.

## Repository layout
```
scripts/ComPlEx/     Comparative co-expression pipeline and conserved-core analysis
scripts/networks/    Per-species whole-genome co-expression networks of the resource
scripts/DE/          Differential expression
scripts/GO/          GO enrichment (Tables S1, S2, S3)
scripts/annotation/  GO-enrichment backgrounds (expressed genes per species)
scripts/figures/     Figure producers
scripts/validation/  Analyses supporting the conserved core
results/             Key reference outputs, for checking a successful run
data/                Not included; download from figshare and unpack here (see Data)
```

`scripts/ComPlEx/` -- the ComPlEx co-expression workflow (Netotea et al. 2014): variance-stabilising
data preparation (`Network_Generation_setup.R`), the 15 pairwise co-expressolog comparisons
(`ComPlEx.R`, run once per species pair; for each pair it builds mutual-rank networks from absolute
Pearson correlation on the genes with an orthologue in that pair and keeps the top 3 % of edges), extraction of the complete cliques and the conserved core (`complete_cliques.R`, `size_rule.R`, `core_cliques.R`, `clique_clusters.R`),
and the co-expressolog table (`Genelists.Rmd`). `submit_ComPlEx_cold.sh` submits `ComPlEx.R` through
`RunR.sh` on a Slurm cluster; adapt both to your environment.

`scripts/networks/build_coexpression_network.R` -- the per-species whole-genome co-expression networks
provided with the resource (`data/networks/network_<species>_hrr100.tsv`), built from raw counts
(`export_network_inputs.R`): genes with at least 1 CPM in at least k samples (k = the smallest replicated
group, at least 3), TMM normalisation, jackknifed Pearson correlation ranked per gene, and an edge kept where
the highest reciprocal rank is at most 100. The method and its R commands are described in
`docs/coexpression_networks.md`. These networks are for exploring a gene's co-expression neighbourhood and are
separate from the pairwise ComPlEx networks.

`scripts/validation/` -- the tests and counts behind the conserved-core results: the label-permutation and
family-size-aware nulls, the STRING cluster-coherence test (Fig. S4) and the counts quoted in the Results. Each is
listed with its inputs and outputs in `REPRODUCE.md`.

## Reproducing the analysis
1. Quantify reads (Salmon) against each species assembly, import with tximport, and VST-normalise
   (DESeq2). Call differential expression (DESeq2) per species and define orthogroups with
   OrthoFinder (as in Rodriguez et al. 2026, *Nat. Commun.* 17, doi:10.1038/s41467-026-75624-2).
2. Compute the 15 pairwise co-expressolog comparisons (`scripts/ComPlEx/ComPlEx.R`), which build
   each pair's mutual-rank networks (top 3 % of edges), and the per-species whole-genome networks of
   the resource (`scripts/networks/export_network_inputs.R`, then `build_coexpression_network.R`).
3. Extract the complete cliques and the conserved core (`complete_cliques.R`, `size_rule.R`, `core_cliques.R`,
   `clique_clusters.R`), and classify each core orthogroup as orthologue- or paralogue-carried (`orthology_class.R`).
4. Produce the figures, tables and numbers of the article from the deposited data: `REPRODUCE.md` lists every
   producer with its inputs and outputs.

Paths in the scripts are given relative to the repository root and expect data under `data/`;
adjust to your local layout where needed.

## Dependencies
R (DESeq2, tximport, tidyverse, here, and the packages listed in `REPRODUCE.md`), with OrthoFinder and Salmon for
the upstream steps.
