#!/usr/bin/env bash
set -euo pipefail
Q="/mnt/c/Users/user/bio-stack/runs/enamel_relics/level4/AMBN_exon3_turtles/queries.fna"
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/level4/AMBN_exon3_turtles"
FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Cpic.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Cpic.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Cmyd.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Cmyd.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Gevg.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Gevg.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Psin.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Psin.tsv" || true
