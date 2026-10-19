#!/usr/bin/env bash
# Passo 5 (WSL, bash -l): busca ML restrita por hipotese e teste AU por fossil.
# Uso: bash -l f04_au.sh [ID ...]   (sem argumento: todas as pastas em au/)
set -euo pipefail
BASE=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis
SEMENTE=20261017
LOG="$BASE/log/comandos_au.txt"
THREADS="${THREADS:-4}"
IDS=("$@")
[ ${#IDS[@]} -eq 0 ] && IDS=($(ls "$BASE/au" | grep -v '\.tsv$'))
for ID in "${IDS[@]}"; do
  M="$BASE/matrizes/$ID"
  for P in "$BASE/au/$ID"/*/; do
    P="${P%/}"
    F="$(basename "$P")"
    : > "$P/ordem.txt"
    : > "$P/hipoteses.trees"
    # arvore ML com esqueleto e todos os fosseis livres (analise 2b)
    REF="$BASE/arvores/$ID/esqueleto.treefile"
    [ "$ID" = MK22 ] && REF="$BASE/arvores/$ID/livre.treefile"
    cat "$REF" >> "$P/hipoteses.trees"
    echo "ML_fosseis_livres" >> "$P/ordem.txt"
    for H in "$P"/H_*.nwk; do
      R="$(basename "$H" .nwk)"
      if [ ! -s "$P/$R.treefile" ]; then
        CMD=(iqtree3 -s "$M/$ID.phy" -st MORPH -p "$M/particoes.nex" -g "$H" --seed "$SEMENTE"
             -T "$THREADS" --prefix "$P/$R" -redo --quiet)
        echo "$(date -Is) ${CMD[*]}" >> "$LOG"
        "${CMD[@]}"
      fi
      cat "$P/$R.treefile" >> "$P/hipoteses.trees"
      echo "${R#H_}" >> "$P/ordem.txt"
    done
    CMD=(iqtree3 -s "$M/$ID.phy" -st MORPH -p "$M/particoes.nex" -te "$REF" -z "$P/hipoteses.trees"
         -n 0 -zb 10000 -zw -au --seed "$SEMENTE" -T "$THREADS" --prefix "$P/AU" -redo --quiet)
    echo "$(date -Is) ${CMD[*]}" >> "$LOG"
    "${CMD[@]}"
    echo "$ID $F: $(wc -l < "$P/ordem.txt") arvores testadas"
  done
done
