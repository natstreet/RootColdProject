#!/bin/bash
# check_start_codons.sh — coordinate check of the promoter windows (promoter_windows.R): the three bases immediately
# downstream of each 1 kb window, on the gene's strand, should be the start codon ATG. Reports per genotype the number of
# windows and how many are followed by ATG. Run in the motif work folder (needs tmp/genomes.tsv, out/<g>_1000.bed,
# samtools on PATH). Usage: bash check_start_codons.sh > out/start_codon_check.tsv
set -euo pipefail
printf "genotype\twindows\tfollowed_by_ATG\n"
while IFS=$'\t' read -r g gff fa pfx rep hard; do
  b=out/${g}_1000.bed
  p=$(awk -F'\t' '$6=="+"{print $1":"$3+1"-"$3+3}' $b | samtools faidx -r /dev/stdin "$fa" | grep -v "^>" | tr a-z A-Z | grep -c "^ATG$" || true)
  m=$(awk -F'\t' '$6=="-"{print $1":"$2-2"-"$2}' $b | samtools faidx -i --mark-strand no -r /dev/stdin "$fa" | grep -v "^>" | tr a-z A-Z | grep -c "^ATG$" || true)
  printf "%s\t%s\t%s\n" "$g" "$(wc -l < $b)" "$((p + m))"
done < tmp/genomes.tsv
