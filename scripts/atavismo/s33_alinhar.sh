#!/usr/bin/env bash
# tBLASTn das consultas reais e embaralhadas do s33 contra o intervalo NRG1..LPL
# de cada especie. Mesmos parametros do s26 (evalue 10, SEG ligado), para que o
# nulo seja o nulo DESTE metodo.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s33_alinhar.sh
set -euo pipefail
A=/mnt/c/Users/Manec/bio-stack/runs/atavismo
O="$A/indice_nulo"
FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'
n=$(grep -c '^>' "$O/consultas.faa")
echo "tblastn $(tblastn -version | head -1 | awk '{print $2}')  |  $n consultas"
for f in "$A"/lncrna/intervalo_*.fna; do
    rot=$(basename "$f" .fna); rot=${rot#intervalo_}
    tblastn -query "$O/consultas.faa" -subject "$f" -evalue 10 -seg yes \
            -max_target_seqs 50 -outfmt "$FMT" > "$O/s33_${rot}.tsv"
    printf '   %-5s %6s HSPs\n' "$rot" "$(wc -l < "$O/s33_${rot}.tsv")"
done
