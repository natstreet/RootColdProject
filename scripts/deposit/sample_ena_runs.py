#!/usr/bin/env python3
"""sample_ena_runs.py — the ENA runs behind every sample, with sequencing depth and platform.

Chain: publication sample code (data/open/samples.tsv) -> sequencing library name (data/quantification/
sample_name_map.tsv) -> the lane file(s) the quantification used (data/quantification/lane_map.tsv, used == yes) ->
the ENA run whose submitted file carries the same library key. The key is the Novogene library prefix
<library>_<EKRN id>_<flowcell>_<lane> for Arabidopsis, aspen and birch, and the facility library number (P#####_###)
for pine and spruce; a spruce library was sequenced as technical replicates on two lanes, and pine sample 120 was
merged from two runs, so those samples have two runs each.

Run metadata come from the ENA portal API (read_count = read pairs, instrument_model, submitted and FASTQ file
names). For each run the first read of each FASTQ file (mate 1 and mate 2) is read from the first 64 kB of the file:
the sequencer's instrument and flowcell fields from its header, and the read length from its sequence. ENA's
instrument_model is the value submitted; the header is what the sequencer wrote. (ENA's base_count is not used: for
some of these runs it counts one mate only.)

Batches: runs are grouped by library naming scheme, and every run in a batch must share one instrument-ID prefix in
its read headers (the script stops otherwise). Batch labels are those of the Methods (sequencing provider and
platform); the platform part is checked against the header instrument prefix in the output, not asserted.

Inputs (deposit v4 or later, unpacked into ./data):
  data/open/samples.tsv, data/quantification/sample_name_map.tsv, data/quantification/lane_map.tsv
Outputs:
  results/ena_run_metadata.tsv   one row per run: run, study, sample_id, library, read_pairs, read_length,
                                 ena_instrument_model, header_instrument, header_flowcell, batch
  results/sample_ena_runs.tsv    one row per sample: sample_id, ena_runs (";"-separated), read_pairs, batch,
                                 runs_of_same_library_not_used (in ENA, but not in the quantification)
  results/sequencing_depth.tsv   per batch: samples, runs, median read pairs per sample, read length(s),
                                 ENA instrument_model value(s), header instrument prefix(es)
Usage: python3 scripts/deposit/sample_ena_runs.py
"""
import csv, gzip, io, os, re, statistics, sys, urllib.parse, urllib.request, zlib
from collections import defaultdict

STUDIES = ["PRJEB104117", "PRJEB104118", "PRJEB104119", "PRJEB104120", "PRJEB26918"]   # Data availability
API = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = "run_accession,study_accession,read_count,instrument_model,submitted_ftp,fastq_ftp"
BATCH_LABEL = {"novogene": "Novogene NovaSeq 6000", "P21452": "NGI NovaSeq 6000", "P1406": "HiSeq 2500"}


def rd(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def ena_runs():
    rows = []
    for st in STUDIES:
        q = urllib.parse.urlencode({"result": "read_run", "query": f"study_accession={st}", "fields": FIELDS,
                                    "format": "tsv", "limit": 0})
        t = urllib.request.urlopen(f"{API}?{q}", timeout=300).read().decode()
        rows += list(csv.DictReader(io.StringIO(t), delimiter="\t"))
    return rows


def first_record(url):
    """header and sequence of the first read in a gzipped FASTQ file, from its first 64 kB"""
    req = urllib.request.Request("https://" + url, headers={"Range": "bytes=0-65535"})
    raw = urllib.request.urlopen(req, timeout=300).read()
    lines = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw).decode(errors="replace").split("\n")
    return lines[0], lines[1]


def sequencer_fields(fastq_ftp):
    """(instrument, flowcell, read lengths) from the sequencer's own header, which ENA keeps as the second
    whitespace-separated field (@<run>.<n> <instrument>:<run>:<flowcell>:<lane>:...)"""
    lens, inst = [], set()
    for url in fastq_ftp.split(";"):
        h, seq = first_record(url)
        f = h.split()[1].split(":")
        inst.add((f[0], f[2])); lens.append(len(seq))
    assert len(inst) == 1, inst
    return inst.pop() + ("+".join(map(str, lens)),)


def key_of(name, library_level=False):
    """run key from a lane-file or submitted-file name: the Novogene prefix up to the lane, or the facility library
    number plus its sample-sheet index (P#####_###_S#); with library_level, the facility library number alone"""
    m = re.search(r"(P\d{4,5}_\d{3})(_S\d+)?", name)
    if m:
        return m.group(1) if library_level or not m.group(2) else m.group(1) + m.group(2)
    m = re.match(r"^(.+?_EKRN\d+-\w+_\w+?_L\d)", name)
    return m.group(1) if m else None


def combines_runs(lane_file):
    """a lane file built from all of its library's runs: a merge, or combined technical replicates"""
    return "_merge_" in lane_file or lane_file.startswith("TechRepComb_")


def batch_of(key):
    if "_EKRN" in key:
        return "novogene"
    return key.split("_")[0]


def main():
    samples = rd("data/open/samples.tsv")
    name = {(r["genotype"], r["deposited_sample_id"]): r["sequencing_name"] for r in rd("data/quantification/sample_name_map.tsv")}
    lanes = defaultdict(list)
    for r in rd("data/quantification/lane_map.tsv"):
        if r["used"] == "yes":
            lanes[(r["genotype"], r["deposited_sample_id"])].append(r["quantification_file"])

    runs = ena_runs()
    by_key, by_lib = defaultdict(list), defaultdict(list)
    for r in runs:
        f = os.path.basename(r["submitted_ftp"].split(";")[0])
        k = key_of(f)
        if k:
            by_key[k].append(r)
            by_lib[key_of(f, library_level=True)].append(r)

    run_rows, sample_rows, used = [], [], set()
    for s in samples:
        sp, sid = s["species"], s["sample_id"]
        lib = name[(sp, sid)]
        keys = {(key_of(f, library_level=combines_runs(f)), combines_runs(f)) for f in lanes[(sp, lib)]}
        if not keys or any(k is None for k, _ in keys):
            sys.exit(f"{sp} {sid}: no library key for lane files {lanes[(sp, lib)]}")
        rr = [r for k, whole in sorted(keys) for r in (by_lib if whole else by_key).get(k, [])]
        libs = {key_of(os.path.basename(r["submitted_ftp"].split(";")[0]), library_level=True) for r in rr}
        left_out = sorted({r["run_accession"] for l in libs for r in by_lib[l]} - {r["run_accession"] for r in rr})
        keys = {k for k, _ in keys}
        if not rr:
            sys.exit(f"{sp} {sid}: no ENA run for keys {sorted(keys)}")
        b = {batch_of(k) for k in keys}
        assert len(b) == 1
        b = b.pop()
        pairs = 0
        for r in rr:
            used.add(r["run_accession"])
            n = int(r["read_count"]); pairs += n
            inst, cell, rlen = sequencer_fields(r["fastq_ftp"])
            run_rows.append({"run": r["run_accession"], "study": r["study_accession"], "sample_id": sid, "species": sp,
                             "library": lib, "read_pairs": n, "read_length": rlen,
                             "ena_instrument_model": r["instrument_model"], "header_instrument": inst,
                             "header_flowcell": cell, "batch": BATCH_LABEL[b]})
        sample_rows.append({"sample_id": sid, "species": sp, "ena_runs": ";".join(sorted(r["run_accession"] for r in rr)),
                            "read_pairs": pairs, "batch": BATCH_LABEL[b],
                            "runs_of_same_library_not_used": ";".join(left_out)})

    # one instrument-ID prefix per batch, from the read headers
    pref = defaultdict(set)
    for r in run_rows:
        pref[r["batch"]].add(re.match(r"[A-Za-z-]+", r["header_instrument"]).group(0))
    for b, p in pref.items():
        if len(p) != 1:
            sys.exit(f"batch {b}: read headers carry more than one instrument prefix {sorted(p)}")

    os.makedirs("results", exist_ok=True)
    for p, rows in (("results/ena_run_metadata.tsv", run_rows), ("results/sample_ena_runs.tsv", sample_rows)):
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    with open("results/sequencing_depth.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["batch", "samples", "runs", "median_read_pairs_per_sample", "read_lengths",
                    "ena_instrument_models", "header_instrument_prefix"])
        for b in sorted({r["batch"] for r in sample_rows}):
            ss = [r for r in sample_rows if r["batch"] == b]; rr = [r for r in run_rows if r["batch"] == b]
            w.writerow([b, len(ss), len(rr), statistics.median(r["read_pairs"] for r in ss),
                        ";".join(sorted({str(r["read_length"]) for r in rr})),
                        ";".join(sorted({r["ena_instrument_model"] for r in rr})), ";".join(sorted(pref[b]))])
    unused = sorted(r["run_accession"] for r in runs if r["run_accession"] not in used)
    print(f"{len(sample_rows)} samples, {len(run_rows)} runs used, {len(unused)} runs in these studies not used: "
          + ", ".join(unused))
    print(open("results/sequencing_depth.tsv").read())


if __name__ == "__main__":
    main()
