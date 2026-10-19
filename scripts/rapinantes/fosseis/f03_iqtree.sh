#!/usr/bin/env bash
# Passo 3 (WSL, bash -l): ML morfologica com IQ-TREE 3, por matriz.
#   livre      morfologia sozinha
#   esqueleto  atuais presos a topologia de Stiller et al. 2024 (-g), fosseis livres
# Uso: bash -l f03_iqtree.sh <ID> [<ID> ...]     (ID = pasta em matrizes/)
# CPU: -T 4 no maximo (outras analises dividem a maquina).
set -euo pipefail
BASE=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis
SEMENTE=20261017
LOG="$BASE/log/comandos_iqtree.txt"
mkdir -p "$BASE/log"
iqtree3 --version | head -1 >> "$LOG"
declare -A EXTERNO=( [K19]=Eudromia_elegans [K17]=Eudromia_elegans [M11]=Tinamidae [M05]=Tinamidae
  [M16]=Anhimidae [MC03]=Apsaravis [LZ07]=Lithornis [VZ26]=Anhimidae [MK22]=Anhimidae [C25]=Gansus_yumenensis [D21]=Meleagris_gallopavo )
for ID in "$@"; do
  M="$BASE/matrizes/$ID"
  O="$BASE/arvores/$ID"
  mkdir -p "$O"
  for MODO in livre esqueleto; do
    EXTRA=()
    if [ "$MODO" = esqueleto ]; then
      [ "$ID" = MK22 ] && continue   # so 3 atuais: esqueleto nao restringe nada
      EXTRA=(-g "$M/esqueleto_stiller.nwk")
    fi
    if [ -s "$O/$MODO.treefile" ] && [ -s "$O/$MODO.ufboot" ]; then
      echo "ja feito: $ID $MODO"; continue
    fi
    CMD=(iqtree3 -s "$M/$ID.phy" -st MORPH -p "$M/particoes.nex" "${EXTRA[@]}"
         -B 1000 --bnni -alrt 1000 -wbt --seed "$SEMENTE" -T "${THREADS:-4}" -o "${EXTERNO[$ID]}"
         --prefix "$O/$MODO" -redo --quiet)
    echo "$(date -Is) ${CMD[*]}" >> "$LOG"
    "${CMD[@]}"
    echo "$ID $MODO: $(grep -m1 'BEST SCORE' "$O/$MODO.log" || true)"
  done
done
