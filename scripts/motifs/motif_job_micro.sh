#!/bin/bash
#SBATCH -J motifs
#SBATCH -p main
#SBATCH --exclude=R2D2,pando,kalkyl
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH -o logs/motifs_%j.out
# Promoter motif track (DECISIONS.md D062 item 4, D063; design R040 as amended in A035), run on micro:
# promoter windows (1 and 2 kb upstream of the ATG) for the genes of motif_pairs.R, sequences and covariates, then FIMO
# (JASPAR 2024 CORE plants non-redundant, p < 1e-4, both strands, 0-order background from each promoter set).
# Genome and annotation per genotype are listed below; the work folder holds scripts/, inputs/ and a pixi environment.
ulimit -c 0
set -euo pipefail
W=/mnt/ada/projects/spruce/nstreet/motifs_20261008
cd $W; mkdir -p logs tmp out
export PIXI_HOME=$W/.pixi_home PIXI_CACHE_DIR=$W/.pixi_cache
if [ ! -x $W/pixi ]; then curl -fsSL https://github.com/prefix-dev/pixi/releases/latest/download/pixi-x86_64-unknown-linux-musl.tar.gz | tar xz -C $W; fi
$W/pixi install
P="$W/pixi run"
$P samtools --version | head -1; $P fimo --version; $P Rscript -e 'cat(R.version.string, "\n")'
R=/mnt/reference
T="/usr/bin/time -f %e_s_%M_kB"
MEME=inputs/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt

# Arabidopsis transposable elements (Araport11) as the repeat intervals. Birch FASTA is gzip, not bgzip, and the spruce
# .fa.gz on micro does not match its .gzi (samtools: "Invalid BGZF header at offset 292948"): both are recompressed here.
gzip -dc $R/Arabidopsis-thaliana/ARAPORT11/gff3/Araport11_GFF3_genes_transposons.201606.gff.gz \
  | awk -F'\t' '$3=="transposable_element" || $3=="transposon_fragment"' > tmp/col0_TE.gff
[ -s tmp/birch.fa.fai ] || { gzip -dc $R/Betula-pendula/v1.4/fasta/genome.fa.gz > tmp/birch.fa; $P samtools faidx tmp/birch.fa; }
[ -s tmp/spruce.fa.gz.gzi ] || { gzip -dc $R/Picea-abies/v2.0/fasta/Picab02_chromosomes_and_unplaced.fa.gz | $P bgzip -@ 8 > tmp/spruce.fa.gz; $P samtools faidx tmp/spruce.fa.gz; }
cmp <(cut -f1,2 tmp/spruce.fa.gz.fai) <(cut -f1,2 $R/Picea-abies/v2.0/fasta/Picab02_chromosomes_and_unplaced.fa.gz.fai)   # same names and lengths

# (pine: the repeat GFF, not the .bed beside it, which is space-separated with a header)
# genotype | GFF3 | FASTA | FASTA name prefix | repeat intervals | hard-masked FASTA
cat > tmp/genomes.tsv <<TSV
col0	$R/Arabidopsis-thaliana/ARAPORT11/gff3/Araport11_GFF3_genes_transposons.201606.gff.gz	$R/Arabidopsis-thaliana/TAIR10/fasta/genome.fa	-	tmp/col0_TE.gff	-
aspen	$R/Populus-tremula/v2.2/gff/Potra02_genes.gff	$R/Populus-tremula/v2.2/fasta/Potra02_genome.fasta	-	-	$R/Populus-tremula/v2.2/fasta/Potra02_genome_hardmasked.fasta
birch	$R/Betula-pendula/v1.4/gff3/Betula_pendula_subsp_pendula_annot.corrected_filtered.gff3	tmp/birch.fa	lcl|Bpe_	-	-
spruce	$R/Picea-abies/v2.0/gff3/Picab02_230926_at01_all_sorted.gff3.gz	tmp/spruce.fa.gz	-	$R/Picea-abies/v2.0/gff3/Picab02_repeats.gff3.gz	-
pine	$R/Pinus-sylvestris/v1.0/gff3/Pinsy01_240308_at01_all.gff3.gz	$R/Pinus-sylvestris/v1.0/fasta/Pinsy01_chromosomes_and_unplaced.fasta.gz	-	$R/Pinus-sylvestris/v1.0/gff3/Pinsy01_chromosomes_and_unplaced.all_repeats.gff.gz	-
TSV
cp tmp/genomes.tsv out/genomes.tsv
while IFS=$'\t' read -r g gff fa pfx rep hard; do
  fai=$fa.fai
  [ -s out/${g}_windows.tsv ] || { echo "== $g windows"; $T $P Rscript scripts/promoter_windows.R "$gff" inputs/genes_$g.txt "$fai" "$pfx" "$rep" out/$g; }
  for win in 1000 2000; do
    [ -s out/${g}_${win}_stats.tsv ] && continue
    echo "== $g $win sequences"
    if [ "$hard" = "-" ]; then $P python3 scripts/extract_promoters.py out/${g}_$win.bed "$fa" out/${g}_$win
    else $P python3 scripts/extract_promoters.py out/${g}_$win.bed "$fa" out/${g}_$win "$hard"; fi
  done
done < tmp/genomes.tsv

# FIMO, one run per genotype and window, 8 at a time
for g in col0 aspen birch spruce pine; do for win in 1000 2000; do echo "$g $win"; done; done > tmp/fimo_runs.txt
fimo_one() {
  g=$1; win=$2; P="$3"; MEME=$4
  $P fasta-get-markov -m 0 out/${g}_$win.fa out/${g}_$win.bg >/dev/null
  $P fimo --thresh 1e-4 --text --skip-matched-sequence --bfile out/${g}_$win.bg $MEME out/${g}_$win.fa 2> logs/fimo_${g}_$win.err | gzip > tmp/fimo_${g}_$win.tsv.gz
  $P python3 scripts/fimo_summary.py tmp/fimo_${g}_$win.tsv.gz out/hits_${g}_$win.tsv.gz
}
export -f fimo_one
$T xargs -P 8 -L 1 bash -c 'fimo_one $0 $1 "'"$P"'" '"$MEME" < tmp/fimo_runs.txt
(cd out && sha256sum * > ../SHA256SUMS_out.txt)
echo done
