# RootColdProject

Analysis code for:

**Comparative genomics of cold-temperature responses in boreal tree roots**
Norway spruce (*Picea abies*), Scots pine (*Pinus sylvestris*), silver birch (*Betula pendula*),
aspen (*Populus tremula*) and two *Arabidopsis thaliana* ecotypes (Col-0, Ost-0).

This repository contains the scripts needed to reproduce the analyses and figures in the paper.
The underlying data are archived separately (see **Data**); download them into a local `data/`
directory before running the scripts.

## Overview
Fine-root transcriptomes of four boreal trees and two *Arabidopsis* ecotypes were profiled over a
ten-day 5 °C (chilling) treatment. Differential expression is largely lineage-specific, but
comparative co-expression (ComPlEx) recovers a conserved regulatory core (116 cliques /
118 orthogroups) shared across all 15 pairwise species comparisons and spanning more than
300 million years of evolution. Supporting analyses test that this core is significant against a
permutation null, is coherent in an independent co-expression compendium (STRING), replicates in
independent data, and shows conserved promoter cis-motifs.

## Data
The scripts read their inputs from a local `data/` directory. Most of it is not held in this
repository and is obtained as below; the one exception is the small per-species super cluster
membership used by Figure S2 and Table S4, which is included directly under
`data/superclusters/membership/`. Obtain the rest as follows:

- **Processed data** (per-species VST expression matrices, sample metadata, DEG tables,
  transcription-factor lists, ComPlEx outputs and clique definitions, functional annotation,
  the OrthoFinder orthogroup table): SciLifeLab figshare, https://figshare.scilifelab.se/
  (DOI:10.17044/scilifelab.32747586). Download `data.tar.gz` and unpack it; it expands to a
  `data/` directory containing `expression/`, `metadata/`, `DEGs/`, `TranscriptionFactors/`,
  `annotation/`, `orthologs/`, `ComPlEx/` and `superclusters/`.
  The `ComPlEx/` step can be re-run either from the cached per-pair orthologue
  and annotation tables (`ComPlEx/orthologs-*.RData`, `ComPlEx/annotation-*.RData`) or, from
  scratch, from `orthologs/orthologs.parquet` together with `annotation/gene_aliases_20140331.txt`
  and `annotation/Orthogroups_130323_predefined_tree.tsv`. The variance-stabilised
  `expression/` matrices are the starting point for network construction; the per-species
  `DESeqDataSet` objects (`dds/`) are included and are the input to
  `scripts/DE/differential_expression.R`.
- **Raw sequencing reads:** ENA umbrella BioProject **PRJEB104158** (components PRJEB104117
  Arabidopsis, PRJEB104118 birch, PRJEB104119 aspen, PRJEB104120 Scots pine, PRJEB26918 Norway
  spruce roots).
- **Genome assemblies:** *Picea abies* v2.0 and *Pinus sylvestris* v1.0 (Kalman et al. 2025),
  *Populus tremula* v2.2, *Betula pendula* v1.4, *A. thaliana* TAIR10 / AtRTD3.

## Repository layout
```
scripts/ComPlEx/     Comparative co-expression pipeline and conserved-core analysis
scripts/validation/  Analyses supporting the conserved core
results/             Key reference outputs, for checking a successful run
data/                Not included; download from figshare and unpack here (see Data)
```

`scripts/ComPlEx/` -- the ComPlEx co-expression workflow (Netotea et al. 2014): variance-stabilising
data preparation (`Network_Generation_setup.R`), per-species mutual-rank network construction
(`ComPlEx_MR_Network.R`), the 15 pairwise co-expressolog comparisons (`ComPlEx.R`, run once per
species pair), extraction of cliques and the conserved core (`cliques_step1.R`, `cliques_step2.R`),
and the downstream visualisation and clique analyses (`1_Data_visualisation.Rmd`,
`2_Sample_comparison.Rmd`, `3_Cliques.Rmd`, `5_Clique_analysis.Rmd`, `Genelists.Rmd`,
`UpsetPlotsEtc.Rmd`). `RunR.sh` / `submitR.sh` / `submit_ComPlEx_cold.sh` are HPC
wrappers; adapt them to your environment.

`scripts/validation/`:
- `conservation_null_test.py` -- permutation null test for the conserved core.
- `string_external_validation.py` -- within-cluster co-expression test in STRING v12.
- `independent_replication.py` -- replication of the core in independent data.
- `cis_motif_conservation.py`, `run_spruce_pine_cismotif.py` -- cross-species promoter cis-motif
  conservation test.

## Reproducing the analysis
1. Quantify reads (Salmon) against each species assembly, import with tximport, and VST-normalise
   (DESeq2). Call differential expression (DESeq2) per species and define orthogroups with
   OrthoFinder (as in Rodriguez et al. 2025).
2. Build per-species mutual-rank co-expression networks (top 3 % of edges) and compute the 15
   pairwise co-expressolog comparisons (`scripts/ComPlEx/`).
3. Extract cliques and the conserved core (`cliques_step1.R`, `cliques_step2.R`).
4. Produce the clique and transcription-factor expression analyses and figures
   (`5_Clique_analysis.Rmd`, `UpsetPlotsEtc.Rmd`, `Genelists.Rmd`).
5. Run the supporting analyses in `scripts/validation/`.

Paths in the scripts are given relative to the repository root and expect data under `data/`;
adjust to your local layout where needed.

## Dependencies
R (DESeq2, tximport, tidyverse, here) and Python 3.10 (numpy, pandas, scipy, statsmodels,
networkx, pyreadr, openpyxl), with OrthoFinder, Salmon and the MEME suite (FIMO) for the upstream
and cis-motif steps.
