#!/bin/bash
# Pontua topologias alternativas fixas (focais, 42 taxa) com ASTRAL-IV (-C) e
# grava o escore total de quartetos de cada uma: a alternativa de maior escore e
# a que o ASTRAL preferiria. Serve para conferir se o q de um ramo (frequencia
# normalizada por gene) e o objetivo (contagem bruta de quartetos) apontam para o mesmo lado.
#   PADRAO='T[78]*' bash -l pontuar_alternativas.sh <tipo> [...]   (PADRAO: glob das topologias, padrao T*)
set -euo pipefail
O=~/rapinantes/reanalise/b10k/aster
T=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/reanalise/b10k/topologias
for tipo in "$@"; do
  for f in "$T"/${PADRAO:-T*}.tre; do
    b=$(basename "$f" .tre)
    astral -C -c "$f" -t 2 --seed 12345 -u 1 -i "$O/$tipo.c95.focal.tre" -o "$O/alt_${tipo}_$b.tre" 2> "$O/alt_${tipo}_$b.log"
    printf '%s\t%s\t%s\n' "$tipo" "$b" "$(grep -E '^Score:' "$O/alt_${tipo}_$b.log" | tail -1 | awk '{print $2}')"
  done
done
