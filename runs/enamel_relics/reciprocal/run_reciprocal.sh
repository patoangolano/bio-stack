#!/usr/bin/env bash
set -euo pipefail
D=/mnt/c/Users/user/bio-stack/runs/enamel_relics
FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle"
blastx -query "$D/reciprocal/candidates.fna" -db "$D/blastdb/amis_proteome" -evalue 10 -seg yes -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$D/reciprocal/blastx_proteome_seg_yes.tsv"
blastx -query "$D/reciprocal/candidates.fna" -db "$D/blastdb/amis_proteome" -evalue 10 -seg no -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$D/reciprocal/blastx_proteome_seg_no.tsv"
blastx -query "$D/reciprocal/candidates_top_hsp.fna" -db "$D/blastdb/amis_proteome" -evalue 10 -seg yes -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$D/reciprocal/blastx_top_hsp_seg_yes.tsv"
# Single-copy test for Gallus gallus psi-ENAM and psi-AMBN against Ggal NRG1-LPL interval
awk '/^>Ggal__ENAM_Amis$|^>Ggal__AMBN_Amis$/{p=1; print; next} /^>/{p=0} p' "$D/reciprocal/candidates.fna" > "$D/reciprocal/ggal_relics.fna"
blastn -task dc-megablast -query "$D/reciprocal/ggal_relics.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-3 -dust yes -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" > "$D/reciprocal/ggal_single_copy.tsv"
