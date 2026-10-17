#!/usr/bin/env bash
# tBLASTn do IBSP do Anolis (real + embaralhadas) contra os intervalos do DSPP do
# s27, que tem o IBSP como ancora. Mesmos parametros do s26/s33.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s34_alinhar.sh
set -euo pipefail
A=/mnt/c/Users/Manec/bio-stack/runs/atavismo
O="$A/composicao"
FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'
for f in "$A"/painel/intervalo_DSPP_*.fna; do
    rot=$(basename "$f" .fna); rot=${rot#intervalo_DSPP_}
    tblastn -query "$O/consultas.faa" -subject "$f" -evalue 10 -seg yes \
            -max_target_seqs 50 -outfmt "$FMT" > "$O/s34_${rot}.tsv"
    printf '   %-5s %5s HSPs\n' "$rot" "$(wc -l < "$O/s34_${rot}.tsv")"
done
