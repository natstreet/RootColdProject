import gzip, re, sys
from collections import defaultdict
import pandas as pd, numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

SD="/sessions/sharp-modest-wrightdata/spruce_pine_sd/"
CATS=SD+"../outputs/complex_out/results/ComPlEx/RData/orthogroup_categories.tsv"
GP=SD+"../outputs/complex_out/results/ComPlEx/RData/weighted_gene_pairs.tsv"
WIN=2000  # proximal-promoter window (bp closest to TSS = 3' end of stored seq)

MOTIFS={"DRE_CRT(CBF)":"RCCGAC","ABRE":"ACGTGKC","EE(circadian)":"AAAATATCT",
        "CBS(CCA1)":"AAGATATTT","G-box":"CACGTG","LTRE":"CCGAAA"}
IU={"A":"A","C":"C","G":"G","T":"T","R":"[AG]","Y":"[CT]","S":"[GC]","W":"[AT]",
    "K":"[GT]","M":"[AC]","N":"[ACGT]"}
rgx={m:re.compile("".join(IU[b] for b in s)) for m,s in MOTIFS.items()}
def rc(s): return s.translate(str.maketrans("ACGTRYKMWSN","TGCAYRMKWSN"))[::-1]
def hit(seq,r): return bool(r.search(seq) or r.search(rc(seq)))

# 1) conserved OGs
cat=pd.read_csv(CATS,sep="\t")
conserved_og=set(cat.loc[cat.category=="conserved","OrthoGroup"])
print(f"conserved OGs: {len(conserved_og)}")

# 2) map genes -> conserved/background per species (PA=spruce col2, PS=pine col3)
gp=pd.read_csv(GP,sep="\t",usecols=["OrthoGroup","Species1","Species2"])
con={"PA":set(),"PS":set()}; allg={"PA":set(),"PS":set()}
for og,a,b in gp.itertuples(index=False):
    allg["PA"].add(a); allg["PS"].add(b)
    if og in conserved_og: con["PA"].add(a); con["PS"].add(b)
bg={sp:allg[sp]-con[sp] for sp in("PA","PS")}
print("spruce conserved/bg genes:",len(con["PA"]),len(bg["PA"]))
print("pine   conserved/bg genes:",len(con["PS"]),len(bg["PS"]))

# 3) load proximal promoters (last WIN bp) for needed genes only
def load_prom(path,needed):
    out={}; name=None; buf=[]
    with gzip.open(path,"rt") as fh:
        for line in fh:
            if line[0]==">":
                if name in needed and buf: out[name]="".join(buf)[-WIN:]
                name=line[1:].strip(); buf=[]
            elif name in needed: buf.append(line.strip())
    if name in needed and buf: out[name]="".join(buf)[-WIN:]
    return out
need_pa=con["PA"]|bg["PA"]; need_ps=con["PS"]|bg["PS"]
prom={"PA":load_prom(SD+"spruce_upstream_20kb.fa.gz",need_pa),
      "PS":load_prom(SD+"pine_upstream_20kb.fa.gz",need_ps)}
print("promoters loaded PA/PS:",len(prom["PA"]),len(prom["PS"]))

# 4) enrichment per species
rows=[]
for sp,lab in (("PA","spruce"),("PS","pine")):
    cg=[g for g in con[sp] if g in prom[sp]]; bgg=[g for g in bg[sp] if g in prom[sp]]
    pv=[]; ms=[]
    for m,r in rgx.items():
        cc=sum(hit(prom[sp][g],r) for g in cg); bc=sum(hit(prom[sp][g],r) for g in bgg)
        OR,p=fisher_exact([[cc,len(cg)-cc],[bc,len(bgg)-bc]],alternative="greater")
        pv.append(p); ms.append(m)
        rows.append(dict(species=lab,motif=m,cons_frac=cc/len(cg),bg_frac=bc/len(bgg),OR=OR,p=p))
    q=multipletests(pv,method="fdr_bh")[1]
    for m,qq in zip(ms,q):
        for d in rows:
            if d["species"]==lab and d["motif"]==m: d["q"]=qq
df=pd.DataFrame(rows)
df.to_csv("sp_cismotif_results.csv",index=False)
print("\n=== cis-motif enrichment (conserved vs background promoters, proximal 2 kb) ===")
print(df[["species","motif","cons_frac","bg_frac","OR","q"]].to_string(index=False,
      float_format=lambda x:f"{x:.3g}"))
enr=df[df.q<0.05]
both=set(enr[enr.species=="spruce"].motif)&set(enr[enr.species=="pine"].motif)
print("\nMotifs enriched in BOTH species (conserved cis-logic, q<0.05):", sorted(both) or "none")
