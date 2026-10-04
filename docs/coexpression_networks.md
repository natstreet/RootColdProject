# Building per-species co-expression networks from RNA-seq counts

This note describes how the per-species co-expression networks hosted in PlantGenIE are built: which genes are kept,
how the counts are normalised, how gene pairs are scored and how the network is cut. Each step is given with the R
code that performs it. The full script is `scripts/networks/build_coexpression_network.R`.

## Input

The method works on one experiment at a time. It needs two tables:

- a count table: raw (not normalised) read counts with genes in rows and samples in columns, gene identifiers in the
  first column;
- a sample table: one row per sample with its experimental group (for example a treatment and time point, or a
  tissue zone).

```r
library(data.table); library(edgeR); library(matrixStats)

t <- fread("counts.tsv", header = TRUE)
g <- t[[1]]; t[[1]] <- NULL
cnt <- rowsum(as.matrix(t), g)                 # rows with the same gene id are summed (transcript-level input)
storage.mode(cnt) <- "double"
groups <- fread("groups.tsv", colClasses = "character")      # columns: sample, group
grp <- setNames(groups$group, groups$sample)
cnt <- cnt[, intersect(colnames(cnt), names(grp))]
cnt <- cnt[, colSums(cnt) > 0]                  # samples without reads are dropped
grp <- grp[colnames(cnt)]
```

## Step 1. Gene filter

A gene is kept when it reaches at least one count per million (CPM) in at least k samples. The number k is tied to
the experimental design: it is the size of the smallest group that has at least two samples, and never less than
three.

```r
group_sizes <- table(grp)
k <- max(3, min(group_sizes[group_sizes >= 2]))
keep <- rowSums(cpm(cnt) >= 1) >= k
x <- cnt[keep, ]
```

The filter removes genes that are silent, or expressed in fewer samples than a single experimental group holds.
Without it, such genes dominate the network. A gene with reads in one sample only correlates perfectly with any other
gene that has reads in the same sample, so these genes form tight clusters with many links. They then appear as the
most connected genes, although their pattern reflects one sample rather than a co-ordinated response. Tying k to the
smallest group keeps genes that are expressed in only one condition, provided that condition is replicated.

## Step 2. Normalisation

Counts are scaled with trimmed mean of M-values (TMM) factors (Robinson and Oshlack, 2010), computed on the kept
genes, and log-transformed.

```r
f <- normLibSizes(DGEList(x), method = "TMM")$samples
eff <- f$lib.size * f$norm.factors                     # effective library sizes
x <- log2(t(t(x) / eff * mean(eff)) + 1)               # counts scaled to the mean library size, then log2(x + 1)
```

Size-factor-adjusted counts gave the most accurate co-expression networks in a comparison of 36 RNA-seq workflows
(Johnson and Krishnan, 2022).

## Step 3. Jackknifed correlation

Each gene pair is scored by a jackknifed Pearson correlation (Heyer et al., 1999). The correlation is recomputed with
each sample left out in turn, and the value closest to zero is kept. When the leave-one-out values differ in sign the
score is zero. A single sample therefore cannot create a strong link on its own: removing an outlying sample removes
the correlation it drives.

The direct computation, for a small number of genes:

```r
jackknife_r <- function(a, b) {
  r <- vapply(seq_along(a), function(s) cor(a[-s], b[-s]), numeric(1))
  if (all(r > 0)) min(r) else if (all(r < 0)) max(r) else 0
}
```

For a whole genome the script computes the same values exactly, but from running sums in blocks of 500 genes. On
centred values, leaving sample s out changes each sum of squares and cross-product by S/(S − 1) times that sample's
term, where S is the number of samples, so no correlation is recomputed from scratch:

```r
X <- x - rowMeans(x); S <- ncol(X); cf <- S / (S - 1)
Cxx <- rowSums(X^2)                                    # sums of squares
b <- 1:500                                             # one block of genes
Cxy <- tcrossprod(X[b, ], X)                           # cross-products of the block with every gene
lo <- matrix(Inf, length(b), nrow(X)); hi <- -lo
for (s in seq_len(S)) {
  vb <- Cxx[b] - cf * X[b, s]^2; va <- Cxx - cf * X[, s]^2
  r <- (Cxy - cf * outer(X[b, s], X[, s])) / sqrt(outer(vb, va))
  r[!is.finite(r)] <- 0                                # a gene constant without sample s
  lo <- pmin(lo, r); hi <- pmax(hi, r)
}
jk <- round(ifelse(lo > 0, lo, ifelse(hi < 0, hi, 0)), 10)
```

Values are rounded to ten decimal places so that values equal up to floating-point error count as ties.

## Step 4. Ranks and highest reciprocal rank

Within each gene's row, its partners are ranked by the jackknifed correlation, strongest positive first (rank 1).
Ties take the average rank and the gene itself is left out. Negative correlations rank last, so the network describes
positive co-expression.

```r
jk[cbind(seq_along(b), b)] <- NA                       # the gene itself
R <- rowRanks(-jk, ties.method = "average")            # 1 = strongest partner
```

Two genes are linked when each is among the other's top partners. The highest reciprocal rank (HRR; Mutwil et al.,
2010) is the larger of the two ranks, HRR(i, j) = max(R[i, j], R[j, i]), and a pair is kept when HRR ≤ 100. Because
both ranks must be at most 100, only each gene's top 100 partners need to be stored. Each kept edge also carries the
mutual rank, sqrt(R[i, j] × R[j, i]) (Obayashi and Kinoshita, 2009).

```r
top <- which(R <= 100, arr.ind = TRUE)                 # per block: each gene's partners within rank 100
d <- data.table(i = b[top[, 1]], j = top[, 2], rank = R[top])
# after all blocks: keep pairs found in both directions
e <- merge(d[i < j], d[i > j, .(i = j, j = i, rank_back = rank)], by = c("i", "j"))
e[, hrr := pmax(rank, rank_back)][, mr := sqrt(rank * rank_back)]
```

HRR cut-offs below 100 give a sparser network, and this can be applied afterwards without rebuilding. In exNet the
density slider does this.

## Output

One row per edge:

| column | meaning |
|---|---|
| source, target | the two genes |
| rank_source_to_target, rank_target_to_source | each gene's rank of the other (1 = strongest) |
| hrr | highest reciprocal rank, the larger of the two ranks |
| mr | mutual rank, the geometric mean of the two ranks |
| jackknife_r | the jackknifed correlation used for ranking |
| pearson_r | the ordinary Pearson correlation of the same normalised values, for reference |

A summary line records the number of samples and groups, k, the genes before and after the filter, the number of
edges and the node degree distribution.

## Why this method

Each choice was compared in 20 experiments (root cold, needle and root cold and drought, and wood development; 1,234
networks). Gene function (GO), stability under resampling of samples, and the number of noise genes among the most
connected genes were scored. A noise gene was defined as expressed in fewer samples than the smallest group, or having
more than half of its variance in one sample. The main results were as follows.

- The gene filter matters most. Without it (keeping every gene with any variation, as in the earlier PlantGenIE
  networks), all 30 of the most connected genes were noise genes in 13 of the 20 experiments. In the median experiment,
  46% of the strongest links (both ranks ≤ 10) involved a noise gene. With the filter the median number of noise genes
  among the 30 hubs was 0, and about 4% of the strongest links involved one.
- The jackknife removed the remaining noise hubs in 8 of the 9 experiments that still had any. It also halved the
  share of strongest links involving a noise gene, to about 2%, at a small cost in stability (about 0.01 on a 0-1
  scale).
- Pearson and Spearman correlation, the rank transformations (HRR, mutual rank, context likelihood) and TMM against
  upper-quartile normalisation differed little in GO recovery. Pearson correlation gave more stable networks than
  Spearman.

## Practical points

- The jackknife multiplies the computation by the number of samples. The 20 PlantGenIE networks (16,000 to 31,000
  genes after filtering, 15 to 105 samples) each took 6 to 54 minutes on one core and 1.3 to 2.6 GB of memory.
- Because both ranks must be at most 100, a gene has at most about 100 links. Degree therefore saturates for the most
  connected genes, and hubs are better compared by the sum or mean of their HRR values.
- The group table should describe the experimental design, because it sets k. If an experiment has no replicated
  groups, k is three.
- Edges with jackknife_r of zero can occur for genes with few strong partners. They rank below every positive partner.
- Networks from different experiments are comparable in their ranks (HRR, mutual rank), not in their correlation
  values.

## Usage

```sh
Rscript scripts/networks/build_coexpression_network.R counts.tsv groups.tsv network_name [HRR_MAX=100] [exclude=s1;s2]
```

The script needs R with the data.table, edgeR and matrixStats packages.

## References

Heyer, L. J., Kruglyak, S., & Yooseph, S. (1999). Exploring expression data: identification and analysis of
coexpressed genes. Genome Research, 9(11), 1106-1115. https://doi.org/10.1101/gr.9.11.1106

Johnson, K. A., & Krishnan, A. (2022). Robust normalization and transformation techniques for constructing gene
coexpression networks from RNA-seq data. Genome Biology, 23, 1. https://doi.org/10.1186/s13059-021-02568-9

Mutwil, M., Usadel, B., Schütte, M., Loraine, A., Ebenhöh, O., & Persson, S. (2010). Assembly of an interactive
correlation network for the Arabidopsis genome using a novel heuristic clustering algorithm. Plant Physiology,
152(1), 29-43. https://doi.org/10.1104/pp.109.145318

Obayashi, T., & Kinoshita, K. (2009). Rank of correlation coefficient as a comparable measure for biological
significance of gene coexpression. DNA Research, 16(5), 249-260. https://doi.org/10.1093/dnares/dsp016

Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq
data. Genome Biology, 11(3), R25. https://doi.org/10.1186/gb-2010-11-3-r25
