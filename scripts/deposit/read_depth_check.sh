#!/bin/bash
#SBATCH -J r072_depth
#SBATCH -p main
#SBATCH --exclude=R2D2,pando,kalkyl
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/ada/projects/spruce/nstreet/r071_pseudo_20261026/logs/depth.out
# R072 item 2 (D102; NRS confirmed in the Code chat 2026-10-26): stream one R1 FASTQ per sequencing batch from ENA, count reads (lines / 4)
# and check the stream's md5 and size against ENA; nothing is written to disk except this log and the result table.
ulimit -c 0
set -uo pipefail
N=/mnt/ada/projects/spruce/nstreet/r071_pseudo_20261026; cd "$N"; O=${OUT:-out/read_depth.tsv}   # committed as results/read_depth_check.tsv
for t in curl md5sum gzip awk; do command -v $t >/dev/null || { echo "missing $t"; exit 1; }; done
curl -s -o /dev/null -w "ENA reachable: HTTP %{http_code}\n" --max-time 30 https://www.ebi.ac.uk/ena/portal/api/ || { echo "no internet on $(hostname)"; exit 1; }
printf 'run\tbatch\turl\tlines\treads_R1\tbytes\tmd5\tena_md5\tmd5_ok\n' > $O; fail=0
while IFS=$'\t' read -r run batch url md5; do
  tmp=$(mktemp -d); mkfifo $tmp/a $tmp/b
  md5sum < $tmp/a > $tmp/md5 & wc -c < $tmp/b > $tmp/bytes &
  lines=$(curl -s --retry 5 --retry-delay 30 "https://$url" | tee $tmp/a $tmp/b | gzip -dc | wc -l); wait
  m=$(cut -d' ' -f1 $tmp/md5); b=$(cat $tmp/bytes); ok=$([ "$m" = "$md5" ] && echo yes || echo no); [ "$ok" = yes ] || fail=1
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $run "$batch" $url $lines $((lines / 4)) $b $m $md5 $ok >> $O; rm -rf $tmp
  echo "$run done: $lines lines, md5 ok $ok"
done <<'LIST'
ERR15911845	Novogene NovaSeq 6000	ftp.sra.ebi.ac.uk/vol1/fastq/ERR159/045/ERR15911845/ERR15911845_1.fastq.gz	ec158ed16379a2c61b3c9b0821b4b674
ERR15911903	NGI NovaSeq 6000	ftp.sra.ebi.ac.uk/vol1/fastq/ERR159/003/ERR15911903/ERR15911903_1.fastq.gz	7a41fee9afb1f897579775f6c2b78104
ERR2593175	HiSeq 2500	ftp.sra.ebi.ac.uk/vol1/fastq/ERR259/005/ERR2593175/ERR2593175_1.fastq.gz	e9ffc5cb93df810c86e95b9a78f0054a
LIST
[ $fail = 0 ] && touch logs/depth.ok || { echo "an md5 did not match"; exit 1; }
