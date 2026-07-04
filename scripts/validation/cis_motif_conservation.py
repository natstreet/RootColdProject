#!/usr/bin/env python3
"""
cis_motif_conservation.py

Cross-species cis-regulatory motif conservation test for the conserved core.

RATIONALE
---------
If the conserved co-expressologs reflect genuinely conserved *regulation* (and
not incidental correlation), their promoters should be enriched for the same
cis-regulatory motifs across species. This is the logic Movahedi et al. (2011)
used to argue that conserved co-expression marks conserved regulatory
interactions. A positive result here is orthogonal, in-silico support for the
"conserved regulatory backbone" claim without any transgenic work.

WHAT IT TESTS
-------------
For each species, and for a panel of known cold/stress/circadian motifs:
  1. Extract promoter sequences (e.g. 1 kb upstream of the TSS) for
     (a) the conserved-core genes and (b) a background gene set.
  2. Scan promoters for each motif (FIMO if available, else a built-in PWM/IUPAC
     scanner).
  3. Test enrichment of each motif in core vs background (Fisher's exact, BH).
  4. Summarise which motifs are enriched in WHICH species, and flag motifs
     enriched across multiple species (= conserved cis-logic).

MOTIF PANEL (edit MOTIFS)
-------------------------
  DRE/CRT  : RCCGAC      (CBF/DREB1 - the core cold element)
  ABRE     : ACGTG[G/T]C (ABA-responsive)
  EE       : AAAATATCT   (Evening Element - circadian)
  CBS      : AAGATATTT   (CCA1-binding site - circadian)
  G-box    : CACGTG
  LTRE     : CCGAAA      (low-temperature-responsive element)
Add/replace with PWMs (MEME format) via MOTIF_MEME_FILE for FIMO scanning.

INPUTS (edit CONFIG)
--------------------
promoters/<species>.fa : promoter FASTA, headers = gene IDs matching the
                         orthogroup/network gene IDs. For spruce and pine these
                         should be extractable from the spruce_pine_sd resource;
                         the other genomes (Arabidopsis, aspen, birch) need
                         promoter sets sourced from their assemblies/GFFs.
core_genes.tsv         : species <tab> gene_id     (conserved-core members)
background_genes.tsv   : species <tab> gene_id     (e.g. all expressed genes)

DEPS: biopython; optional MEME-suite `fimo` on PATH for PWM scanning.
Run:  python cis_motif_conservation.py
"""
from __future__ import annotations
import os, glob, subprocess, shutil, re
from collections import defaultdict
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ------------- CONFIG -------------
PROMOTER_DIR     = "promoters"
CORE_FILE        = "core_genes.tsv"
BACKGROUND_FILE  = "background_genes.tsv"
MOTIF_MEME_FILE  = None        # set to a MEME-format PWM file to use FIMO instead of IUPAC
FIMO_QVAL        = 1e-3
MIN_SPECIES_FOR_CONSERVED = 3  # motif enriched in >= this many species = "conserved"
OUT_PREFIX       = "cis_motif"

MOTIFS = {                     # IUPAC consensus (used if MOTIF_MEME_FILE is None)
    "DRE_CRT": "RCCGAC",
    "ABRE":    "ACGTGKC",
    "EE":      "AAAATATCT",
    "CBS":     "AAGATATTT",
    "Gbox":    "CACGTG",
    "LTRE":    "CCGAAA",
}
# ----------------------------------

IUPAC = {"A":"A","C":"C","G":"G","T":"T","R":"[AG]","Y":"[CT]","S":"[GC]",
         "W":"[AT]","K":"[GT]","M":"[AC]","B":"[CGT]","D":"[AGT]","H":"[ACT]",
         "V":"[ACG]","N":"[ACGT]"}

def iupac_regex(motif): return "".join(IUPAC[b] for b in motif.upper())

def revcomp(seq):
    return seq.translate(str.maketrans("ACGTRYSWKMBDHVN","TGCAYRSWMKVHDBN"))[::-1]

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]; buf = []
            else:
                buf.append(line.strip().upper())
    if name: seqs[name] = "".join(buf)
    return seqs

def has_motif_iupac(seq, rgx):
    return bool(re.search(rgx, seq) or re.search(rgx, revcomp(seq)))

def load_sets(path):
    d = defaultdict(set)
    for r in pd.read_csv(path, sep="\t", names=["species","gene"], header=0).itertuples():
        d[r.species].add(r.gene)
    return d

def main():
    core = load_sets(CORE_FILE)
    bg   = load_sets(BACKGROUND_FILE)
    use_fimo = MOTIF_MEME_FILE and shutil.which("fimo")
    print(f"Scanning mode: {'FIMO/PWM' if use_fimo else 'built-in IUPAC'}")

    rows = []
    for fa in sorted(glob.glob(os.path.join(PROMOTER_DIR, "*.fa"))):
        sp = os.path.splitext(os.path.basename(fa))[0]
        proms = read_fasta(fa)
        core_g = [g for g in core.get(sp, []) if g in proms]
        bg_g   = [g for g in bg.get(sp, [])   if g in proms and g not in set(core_g)]
        if not core_g or not bg_g:
            print(f"  {sp}: insufficient promoters (core={len(core_g)}, bg={len(bg_g)}) - skipped")
            continue
        rgx = {m: iupac_regex(c) for m, c in MOTIFS.items()}  # IUPAC mode
        pvals, names = [], []
        for m in MOTIFS:
            cc = sum(has_motif_iupac(proms[g], rgx[m]) for g in core_g)
            bc = sum(has_motif_iupac(proms[g], rgx[m]) for g in bg_g)
            table = [[cc, len(core_g)-cc], [bc, len(bg_g)-bc]]
            OR, p = fisher_exact(table, alternative="greater")
            pvals.append(p); names.append(m)
            rows.append({"species": sp, "motif": m, "core_frac": cc/len(core_g),
                         "bg_frac": bc/len(bg_g), "odds_ratio": OR, "p": p})
        rej, q, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
        for nm, qq, rj in zip(names, q, rej):
            for row in rows:
                if row["species"] == sp and row["motif"] == nm:
                    row["q"] = qq; row["enriched"] = bool(rj)
        print(f"  {sp}: {len(core_g)} core / {len(bg_g)} bg promoters scanned")

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT_PREFIX}_per_species.csv", index=False)

    # conservation summary: motifs enriched across species
    enr = df[df.get("enriched", False) == True]
    summary = (enr.groupby("motif")["species"].nunique()
                  .reset_index(name="n_species_enriched")
                  .sort_values("n_species_enriched", ascending=False))
    summary["conserved"] = summary["n_species_enriched"] >= MIN_SPECIES_FOR_CONSERVED
    summary.to_csv(f"{OUT_PREFIX}_conservation_summary.csv", index=False)
    print("\nConserved cis-logic (motifs enriched in >= "
          f"{MIN_SPECIES_FOR_CONSERVED} species):")
    print(summary.to_string(index=False))
    print(f"\nWrote {OUT_PREFIX}_per_species.csv and {OUT_PREFIX}_conservation_summary.csv")

if __name__ == "__main__":
    main()
