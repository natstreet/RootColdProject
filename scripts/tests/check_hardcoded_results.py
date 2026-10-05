#!/usr/bin/env python3
"""check_hardcoded_results.py — find result values written as literals in the analysis code.

A result (a count, a statistic) must be read at run time from a producer's output, never typed into a script: a
literal silently goes stale when the data change (FAILURES.md F010: the core size 132 was fixed in
core_module_directions.R and broke when the core became 137).

The values to look for are read from one or more value tables (TSV): every cell that parses as a number, from the
columns named with --columns. Integers below --min-int are skipped (small integers such as 2, 6 or 15 are mostly
design constants: genotypes, comparisons, thresholds). Each code line under scripts/ is searched for those numbers as
whole tokens, with comments removed (text after '#', outside strings, for R, Python and shell). A hit is a lead to
inspect, not proof of a fault: a constant may coincide with a result.

    python3 scripts/tests/check_hardcoded_results.py results/values_published_R028_R030.tsv \
        [--columns published R028 R030] [--min-int 20] [--allow scripts/tests/hardcoded_allow.tsv]

The allow list (TSV: file, value, reason) records constants inspected and found not to be results; a hit is allowed
when its file and value are listed. Exit status 1 when any hit is not allowed.
"""
import argparse, csv, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("tables", nargs="+")
ap.add_argument("--columns", nargs="+", default=["published", "R028", "R030"])
ap.add_argument("--min-int", type=int, default=20)
ap.add_argument("--allow", default="scripts/tests/hardcoded_allow.tsv")
ap.add_argument("--root", default="scripts")
a = ap.parse_args()

values = set()
for t in a.tables:
    with open(t, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for c in a.columns:
                v = (row.get(c) or "").strip()
                if not re.fullmatch(r"-?\d+(\.\d+)?", v):
                    continue
                if "." not in v and abs(int(v)) < a.min_int:
                    continue
                if "." in v and len(v.replace("-", "").replace(".", "").lstrip("0")) < 3:
                    continue                                  # 0.5, 1.2: too generic to flag
                values.add(v)
if not values:
    sys.exit("no values read")

files = subprocess.run(["git", "ls-files", a.root], capture_output=True, text=True, check=True).stdout.split()
files = [f for f in files if f.endswith((".R", ".Rmd", ".py", ".sh")) and "/reference/" not in f]


def code_part(line):
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        elif ch in "\"'":
            q = ch; out.append(ch)
        elif ch == "#":
            break
        else:
            out.append(ch)
    return "".join(out)


pat = re.compile(r"(?<![\w.])(" + "|".join(sorted(map(re.escape, values), key=len, reverse=True)) + r")(?![\w.])")
hits = []
for f in files:
    with open(f, errors="replace") as fh:
        for i, line in enumerate(fh, 1):
            for m in pat.finditer(code_part(line)):
                hits.append((f"{f}:{i}", m.group(1), line.strip()[:140]))
allowed = set()
if a.allow and os.path.exists(a.allow):
    with open(a.allow, newline="") as fh:
        allowed = {(r["file"], r["value"]) for r in csv.DictReader(fh, delimiter="\t")}
ok = lambda loc, v: (loc.rsplit(":", 1)[0], v) in allowed
bad = [h for h in hits if not ok(h[0], h[1])]
for loc, v, text in hits:
    print(f"{'allowed' if ok(loc, v) else 'HIT'}\t{loc}\t{v}\t{text}")
print(f"{len(values)} result values searched in {len(files)} files: {len(hits)} hits, {len(bad)} not allowed")
sys.exit(1 if bad else 0)
