#!/usr/bin/env bash
set -euo pipefail
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/v102/ambn_exon4_coverage_audit"
INT="/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel"
> "$O/blastn_ws7_all.tsv"
> "$O/tblastn_ws2_all.tsv"
for f in "$INT"/*.fna; do
  sp=$(basename "$f" .fna)
  blastn -task blastn -word_size 7 -evalue 10 -dust no -query "$O/q_ex4_nt.fna" -subject "$f" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | awk -v s="$sp" '{print s"\t"$0}' >> "$O/blastn_ws7_all.tsv" || true
  tblastn -word_size 2 -evalue 50 -seg no -comp_based_stats 0 -query "$O/q_ex4_aa.faa" -subject "$f" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | awk -v s="$sp" '{print s"\t"$0}' >> "$O/tblastn_ws2_all.tsv" || true
done
