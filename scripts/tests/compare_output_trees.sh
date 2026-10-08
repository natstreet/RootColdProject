#!/bin/bash
# compare_output_trees.sh — byte comparison of the producer outputs of two working trees (ledger R042 item 2: which outputs
# of the re-run differ from the previous run). Compares files in the top level, results/ and tables/ (not data/, logs/).
# PDFs are also compared with their CreationDate/ModDate lines removed, since those change on every run.
# Output (stdout): status<TAB>file, status = identical | identical_except_pdf_dates | differs | only_old | only_new
# Usage: bash compare_output_trees.sh <old tree> <new tree>
set -euo pipefail
O=$1; N=$2
list() { (cd "$1" && find . -maxdepth 1 -type f; find ./results ./tables -type f 2>/dev/null) | sed 's|^\./||' | sort; }
comm -3 <(list "$O") <(list "$N") | while IFS= read -r l; do case "$l" in $'\t'*) echo -e "only_new\t${l#$'\t'}";; *) echo -e "only_old\t$l";; esac; done
comm -12 <(list "$O") <(list "$N") | while IFS= read -r f; do
  if cmp -s "$O/$f" "$N/$f"; then echo -e "identical\t$f"
  elif [[ "$f" == *.pdf ]] && cmp -s <(LC_ALL=C grep -av "CreationDate\|ModDate" "$O/$f") <(LC_ALL=C grep -av "CreationDate\|ModDate" "$N/$f"); then echo -e "identical_except_pdf_dates\t$f"
  else echo -e "differs\t$f"; fi
done
