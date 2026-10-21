#!/bin/bash
# check_committed_outputs.sh — the repository's committed derived files against the deposit and a fresh producer run
# (ledger F034: tables committed before the rebuild stayed in results/ and data/ after it). Run before every release.
#   1. Every committed file under data/ that the deposit also holds is byte-identical to the deposit's copy.
#   2. Every committed file under results/ is byte-identical to the output of the same name in <run dir>, a working
#      tree where the producers were run on the deposit's data (scripts/reproduce_figures_tables.sh or the producers
#      named in REPRODUCE.md). Producers write either to the working directory or to results/; both are searched.
#      results/ena_run_metadata.tsv is compared row-sorted: its row order follows the ENA portal's response.
#   A committed file with no counterpart fails.
# Usage (repository root): bash scripts/tests/check_committed_outputs.sh <deposit data dir> <run dir>
set -uo pipefail
DEP=${1:?deposit data dir}; RUN=${2:?run dir}
fail=0; n_data=0; n_res=0
for f in $(git ls-files data); do
  [ -f "$DEP/${f#data/}" ] || continue
  n_data=$((n_data + 1))
  cmp -s "$f" "$DEP/${f#data/}" || { echo "DIFFERS from the deposit: $f"; fail=1; }
done
for f in $(git ls-files results); do
  n_res=$((n_res + 1)); rel=${f#results/}; out=""
  for c in "$RUN/$rel" "$RUN/results/$rel" "$RUN/$(basename "$rel")"; do [ -f "$c" ] && { out=$c; break; }; done
  if [ -z "$out" ]; then echo "NO PRODUCER OUTPUT in the run dir: $f"; fail=1; continue; fi
  if [ "$rel" = ena_run_metadata.tsv ]; then
    cmp -s <(sort "$f") <(sort "$out") || { echo "DIFFERS from the fresh run (rows): $f"; fail=1; }
  else
    cmp -s "$f" "$out" || { echo "DIFFERS from the fresh run: $f"; fail=1; }
  fi
done
echo "data/ files also in the deposit: $n_data; results/ files: $n_res"
[ "$fail" -eq 0 ] && echo "committed outputs: no mismatches" || { echo "committed outputs: mismatches above"; exit 1; }
