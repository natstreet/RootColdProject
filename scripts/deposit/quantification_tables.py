#!/usr/bin/env python3
"""quantification_tables.py — the gene-level tables of the deposit's data/quantification/ from its per-run Salmon files.

For each genotype, every sequencing run of a sample is used (the deposited DESeq2 counts use every run; ledger F032 records
that the earlier tables used one of OstR1001's two runs). Writes, per genotype:
  <g>_gene_counts.tsv.gz   per run, transcript reads summed per gene and rounded; runs then added (the convention of the
                           deposited counts, which DESeq2's collapseReplicates also follows)
  <g>_gene_tpm.tsv.gz      runs combined per transcript (reads summed; effective length the read-weighted mean, the plain
                           mean when the transcript has no reads); TPM recomputed over the transcripts of the genes kept,
                           then summed per gene; 4 decimals
  <g>_gene_efflen.tsv.gz   per gene, the read-weighted mean effective length of its transcripts (plain mean without
                           reads); 2 decimals
and lane_map.tsv (every run file and the sample it belongs to, all used).

Genes kept: the rows of data/expression/<g>_root_cold_stress_expression.txt (the genes of the deposited DESeq2 object);
transcripts of other targets (spruce: transposable elements and repeats in the index) are left out and TPM is rescaled
over the rest. Gene id = transcript id without its suffix (Arabidopsis and aspen ".N", birch ".mN", spruce and pine
".mRNA.N"). Samples, and which run file belongs to which, come from the existing lane_map.tsv (its `used` column is
ignored); columns are in the order of that file. Checking these counts against the deposited DESeq2 counts is
check_quantification_tables.R.

    python3 quantification_tables.py <deposit root> <out dir>
"""
import csv, gzip, os, re, sys
from collections import OrderedDict, defaultdict

ROOT, OUT = sys.argv[1:3]
Q = os.path.join(ROOT, "data/quantification")
SUFFIX = {"col0": r"\.[0-9]+$", "ost0": r"\.[0-9]+$", "aspen": r"\.[0-9]+$", "birch": r"\.m[0-9]+$",
          "spruce": r"\.mRNA\.[0-9]+$", "pine": r"\.mRNA\.[0-9]+$"}
os.makedirs(OUT, exist_ok=True)

lanes = list(csv.DictReader(open(os.path.join(Q, "lane_map.tsv"), newline=""), delimiter="\t"))
samples = OrderedDict()
for r in lanes:
    samples.setdefault(r["genotype"], OrderedDict()).setdefault(r["deposited_sample_id"], []).append(r["quantification_file"])


def read_quant(f):
    with gzip.open(os.path.join(Q, "lanes", f + ".quant.sf.gz"), "rt", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t"); h = next(rd)
        i_n, i_e, i_r = h.index("Name"), h.index("EffectiveLength"), h.index("NumReads")
        return [(r[i_n], float(r[i_e]), float(r[i_r])) for r in rd]


def num(x, d):
    return str(round(x, d))          # Python's float text, as the tables were written ("0.0", "30.7906")


lane_rows = []
for g, smp in samples.items():
    with open(os.path.join(ROOT, f"data/expression/{g}_root_cold_stress_expression.txt"), newline="") as fh:
        rd = csv.reader(fh, delimiter="\t"); next(rd)
        genes = [r[0] for r in rd]
    keep = set(genes); sfx = re.compile(SUFFIX[g])
    counts, tpm, efflen = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    left_out = None
    for s, files in smp.items():
        reads, wlen, lens, n_runs = defaultdict(float), defaultdict(float), defaultdict(list), 0
        gene_rounded = defaultdict(int)
        for f in files:
            q = read_quant(f); n_runs += 1
            per_gene = defaultdict(float)
            for name, e, r in q:
                reads[name] += r; wlen[name] += e * r; lens[name].append(e)
                per_gene[sfx.sub("", name)] += r
            for gid, v in per_gene.items():
                gene_rounded[gid] += int(round(v))
            lane_rows.append((g, s, f, "yes", "single run" if len(files) == 1 else f"summed over {len(files)} runs"))
        rate, eff = {}, {}
        for name in reads:
            eff[name] = wlen[name] / reads[name] if reads[name] > 0 else sum(lens[name]) / len(lens[name])
            rate[name] = reads[name] / eff[name] if eff[name] > 0 else 0.0
        kept = [t for t in reads if sfx.sub("", t) in keep]
        if left_out is None:
            left_out = len(reads) - len(kept)
        tot = sum(rate[t] for t in kept)
        gt, gr, gw, gl = defaultdict(float), defaultdict(float), defaultdict(float), defaultdict(list)
        for t in kept:
            gid = sfx.sub("", t)
            gt[gid] += 1e6 * rate[t] / tot; gr[gid] += reads[t]; gw[gid] += eff[t] * reads[t]; gl[gid].append(eff[t])
        for gid in genes:
            counts[gid][s] = gene_rounded.get(gid, 0)
            tpm[gid][s] = gt.get(gid, 0.0)
            efflen[gid][s] = (gw[gid] / gr[gid] if gr.get(gid, 0) > 0 else (sum(gl[gid]) / len(gl[gid]) if gl.get(gid) else 0.0))
    cols = list(smp)
    for kind, tab, fmt in (("counts", counts, lambda v: str(v)), ("tpm", tpm, lambda v: num(v, 4)), ("efflen", efflen, lambda v: num(v, 2))):
        with gzip.open(os.path.join(OUT, f"{g}_gene_{kind}.tsv.gz"), "wt", newline="") as fh:
            fh.write("gene\t" + "\t".join(cols) + "\n")
            for gid in genes:
                fh.write(gid + "\t" + "\t".join(fmt(tab[gid][s]) for s in cols) + "\n")
    print(f"{g}: {len(genes)} genes, {len(cols)} samples, {sum(len(v) for v in smp.values())} run files; transcripts left out per run: {left_out}")

with open(os.path.join(OUT, "lane_map.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genotype", "deposited_sample_id", "quantification_file", "used", "note"])
    w.writerows(lane_rows)
print(f"lane_map.tsv: {len(lane_rows)} run files, all used")
