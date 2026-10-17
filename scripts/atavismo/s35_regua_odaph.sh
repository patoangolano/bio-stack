#!/usr/bin/env bash
# Regua de sensibilidade do s35: com que bitscore a ODAPH do jacare acha uma
# ODAPH FUNCIONAL a distancia maior que a das aves?
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s35_regua_odaph.sh
#
# Anolis e Podarcis (Lepidosauria) estao a ~280 Ma do jacare; as aves, a ~245.
# Se a ODAPH do jacare acha a deles com folga acima do piso de ruido (22-30
# bits; os "hits" das aves no genoma inteiro ficaram em 32-33), um ODAPH real
# nas aves apareceria acima disso — e o negativo do genoma inteiro vale.
set -euo pipefail
O=/mnt/c/Users/Manec/bio-stack/runs/atavismo/odam_odaph
FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'
for rot in Acar Pmur Cpor; do
    j="$O/janela_odaph_${rot}.fna"
    [ -s "$j" ] || { echo "$rot: sem janela"; continue; }
    tblastn -query "$O/modelo_odaph_Amis.faa" -subject "$j" -seg yes -evalue 10 -outfmt "$FMT" |
        sort -t$'\t' -k10,10gr | head -3 |
        awk -F'\t' -v r="$rot" '{printf "%-5s q%s-%s  %5s%% id  bits %6s  e=%s\n", r, $5, $6, $3, $10, $9}'
done
