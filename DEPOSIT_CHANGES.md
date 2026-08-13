# Deposit / repository clarifications (DEG-definition naming)

Written for anyone reusing the FigShare deposit or this repository. **No manuscript number,
figure or data value changed** as a result of any of this. The manuscript's stated DEG definition
(adjusted P ≤ 0.05 AND |log2FC| ≥ 1, a two-fold change) is correct, and the published figures are
correct. What follows is documentation/naming hygiene, because the artefacts as originally
committed could be — and were — misread.

## What was ambiguous

An independent pre-submission reviewer, working only from the deposit and committed scripts,
concluded the DEG definition in the Methods did not describe the analysis that made the figures,
and flagged it as the paper's most serious problem. That conclusion was wrong but reasonable: the
committed `make_Fig2.R` header said the DEGs used "NO fold-change filter", the per-timepoint lists
are named `_lfc0` (reads as "log-fold-change 0" = no filter), and the `_lfc0` list lengths match
the Figure 2 bars — apparently confirming it. In fact the `_lfc0` lists *are* the two-fold DEGs.

## What each artefact actually is (all verified with commands; see check script below)

- **`DEGs/per_timepoint_DEG_lists/<sp>/…_p0.05_{lfc0,lfc0585,l2fc1}.txt`** — the suffix names the
  `lfcThreshold` argument to DESeq2 `results()`, NOT the fold cut. **All three families carry a
  post-hoc |log2FC| ≥ 1 filter**; every gene in every list is two-fold. `_lfc0` (lfcThreshold = 0)
  is the family used throughout the paper; `_lfc0585` and `_l2fc1` are stricter alternatives.
- **`DEGs/DEGs_<species>.RData`** — the per-species union of the `_lfc0` lists (Col-0 5689, Ost-0
  7995, aspen 11729, birch 3312, spruce 6353, pine 4717). These feed `og_summary` and Figures 1–2.
- **`DEGs/padj_<code>.RData`** — adjusted P-values from the **`lfcThreshold = 1`** test.
  Thresholding at padj ≤ 0.05 reproduces the strict `_l2fc1` lists (24/24 species × timepoint),
  NOT the `_lfc0` DEGs. Do not use it to recheck a paper DEG or a Figure 6 significance mark.
- **`DEGs/log2FoldChange_ps_n.RData`** — a second pine log2FC object of undetermined provenance;
  differs substantially from `log2FoldChange_ps.RData` (the published pine object). Not used by the paper.

## What changed here

- Corrected the misleading "NO fold-change filter" wording in `scripts/figures/make_Fig2.R`, and
  aligned the notes in `scripts/DE/differential_expression.R` and `REPRODUCE.md`
  (commit: "Clarify DEG-list naming in producer docs").
- Expanded `scripts/DE/README.md` to state all of the above once, plainly.
- Added `scripts/DE/check_deposit_consistency.R` (asserts the invariants; PASS/FAIL).

## What was deliberately left alone

- **No deposited data file was renamed or altered.** The `_lfc0` naming stays (a rename would break
  anyone who has downloaded or cited the deposit); it is documented instead. `padj_*` and `_ps_n`
  stay in place, documented, pending any future versioned deposit.
- The manuscript is untouched by this task.

## Could not resolve

- The exact provenance of `log2FoldChange_ps_n.RData` (different normalisation? superseded run?)
  could not be established from the deposited files alone. Documented as non-published; recommend
  removal in a future deposit version once its origin is confirmed.
