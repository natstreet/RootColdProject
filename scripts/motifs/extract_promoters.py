#!/usr/bin/env python3
"""extract_promoters.py — promoter sequences and sequence covariates for the promoter motif track (D062 item 4, D063).

Reads a BED6 of promoter windows (promoter_windows.R; name = gene ID), extracts each window from an indexed FASTA with
`samtools faidx` (minus-strand windows reverse-complemented with -i), and writes:
  <out>.fa         sequences in upper case (the input to FIMO)
  <out>_stats.tsv  gene, length, gc_frac (G+C over non-N bases), n_frac, lower_frac (soft-masked share, meaningful only if
                   the FASTA is soft-masked) and, if a second (hard-masked) FASTA is given, hardmask_frac = share of bases
                   that are N in the hard-masked copy but not in the unmasked one
Usage: extract_promoters.py <bed> <fasta> <out prefix> [<hard-masked fasta>]
"""
import subprocess, sys, tempfile, os

bed, fa, out = sys.argv[1:4]
hard = sys.argv[4] if len(sys.argv) > 4 else None
rows = [l.rstrip("\n").split("\t") for l in open(bed) if l.strip()]

def fetch(fasta, rows):
    seqs = {}
    for strand in "+-":
        sub = [r for r in rows if r[5] == strand]
        if not sub:
            continue
        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".txt") as fh:
            for r in sub:
                fh.write(f"{r[0]}:{int(r[1]) + 1}-{r[2]}\n")
            reg = fh.name
        cmd = ["samtools", "faidx", "-r", reg, fasta] + (["-i", "--mark-strand", "no"] if strand == "-" else [])
        txt = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        os.unlink(reg)
        recs = txt.split(">")[1:]
        if len(recs) != len(sub):
            sys.exit(f"samtools returned {len(recs)} records for {len(sub)} {strand} windows")
        for r, rec in zip(sub, recs):
            head, *lines = rec.split("\n")
            if head.strip() != f"{r[0]}:{int(r[1]) + 1}-{r[2]}":
                sys.exit(f"record order mismatch: {head} vs {r[:3]}")
            seqs[r[3]] = "".join(lines)
    return seqs

s = fetch(fa, rows)
h = fetch(hard, rows) if hard else None
with open(out + ".fa", "w") as fo, open(out + "_stats.tsv", "w") as so:
    so.write("gene\tlength\tgc_frac\tn_frac\tlower_frac" + ("\thardmask_frac" if h else "") + "\n")
    for r in rows:
        g = r[3]; x = s[g]; L = len(x)
        if L != int(r[2]) - int(r[1]):
            sys.exit(f"length mismatch for {g}")
        up = x.upper(); n = up.count("N"); acgt = L - n
        gc = (up.count("G") + up.count("C")) / acgt if acgt else float("nan")
        low = sum(c.islower() for c in x) / L if L else float("nan")
        line = f"{g}\t{L}\t{gc:.4f}\t{n / L:.4f}\t{low:.4f}"
        if h:
            y = h[g].upper()
            if len(y) != L:
                sys.exit(f"hard-masked length mismatch for {g}")
            line += f"\t{sum(1 for a, b in zip(up, y) if b == 'N' and a != 'N') / L:.4f}"
        so.write(line + "\n")
        fo.write(f">{g}\n{up}\n")
print(f"{len(rows)} windows -> {out}.fa")
