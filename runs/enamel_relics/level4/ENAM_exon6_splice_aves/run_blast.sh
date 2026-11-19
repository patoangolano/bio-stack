#!/usr/bin/env bash
set -euo pipefail
Q="/mnt/c/Users/user/bio-stack/runs/enamel_relics/level4/ENAM_exon6_splice_aves/queries.fna"
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/level4/ENAM_exon6_splice_aves"
FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Ggal.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Ggal.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Apla.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Apla.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Shab.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Shab.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Tgut.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Tgut.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Dnov.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Dnov.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Scam.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Scam.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Arow.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Arow.tsv" || true
blastn -task blastn -word_size 7 -query "$Q" -subject "/mnt/c/Users/user/bio-stack/runs/enamel_relics/intervals_enamel/interval_Tgua.fna" -evalue 1e-4 -outfmt "$FMT" > "$O/blast_Tgua.tsv" || true
