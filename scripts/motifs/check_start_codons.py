#!/usr/bin/env python3
"""check_start_codons.py — per-gene start-codon flag for the promoter windows (promoter_windows.R; D064).

The first three coding bases of each window's transcript (column first3 of <genotype>_windows.tsv: genomic pieces in
transcript order, so a start codon split by an intron is read correctly) are fetched on the gene's strand with
`samtools faidx` (-i for minus) and joined. Genes whose model has no annotated start codon (not ATG) are excluded from the
primary motif test and included in its sensitivity run (D064). The window itself is anchored on the first coding base.
Run in the motif work folder (tmp/genomes.tsv, out/<genotype>_windows.tsv; samtools on PATH).
Output (stdout): genotype, gene, transcript, strand, codon, start_codon (TRUE if ATG).
Usage: python3 check_start_codons.py > out/start_codons.tsv
"""
import csv, subprocess, sys

def fetch(fa, regions, minus):
    if not regions:
        return []
    cmd = ["samtools", "faidx", fa, *regions] + (["-i", "--mark-strand", "no"] if minus else [])
    txt = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    recs = [r.split("\n", 1) for r in txt.split(">")[1:]]
    if len(recs) != len(regions) or any(h.strip() != r for (h, _), r in zip(recs, regions)):
        sys.exit("samtools records do not match the requested regions")
    return [s.replace("\n", "").upper() for _, s in recs]

print("genotype\tgene\ttranscript\tstrand\tcodon\tstart_codon")
for g, gff, fa, *_ in csv.reader(open("tmp/genomes.tsv"), delimiter="\t"):
    rows = [r for r in csv.DictReader(open(f"out/{g}_windows.tsv"), delimiter="\t") if r["win"] == "1000"]
    for strand in "+-":
        sub = [r for r in rows if r["strand"] == strand]
        pieces = [p for r in sub for p in r["first3"].split(",")]
        seqs = iter(fetch(fa, pieces, strand == "-"))
        for r in sub:
            codon = "".join(next(seqs) for _ in r["first3"].split(","))
            print("\t".join([g, r["gene"], r["transcript"], strand, codon, "TRUE" if codon == "ATG" else "FALSE"]))
