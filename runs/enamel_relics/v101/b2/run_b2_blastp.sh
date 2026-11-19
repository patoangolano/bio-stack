#!/usr/bin/env bash
set -euo pipefail
D="/mnt/c/Users/user/bio-stack/runs/enamel_relics"
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/b2"
blastp -query "$O/b2_xp_proteins.faa" -db "$D/blastdb/amis_proteome" -evalue 10 -seg yes -max_target_seqs 10 -num_threads 8 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" > "$O/b2_blastp_vs_amis.tsv"
blastp -query "$O/b2_xp_proteins.faa" -db "$D/blastdb/amis_proteome" -evalue 10 -seg no -max_target_seqs 10 -num_threads 8 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" > "$O/b2_blastp_vs_amis_segno.tsv"
