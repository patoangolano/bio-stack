#!/usr/bin/env bash
# AU so entre variantes do esqueleto molecular (sem a arvore morfologica livre,
# que domina o conjunto e rejeita todas): a morfologia prefere o Falconidae
# junto dos gavioes a posicao molecular dele?
# Uso: bash -l f05_au_esqueletos.sh ID [ID ...]
set -euo pipefail
BASE=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis
SEMENTE=20261017
LOG="$BASE/log/comandos_au.txt"
for ID in "$@"; do
  M="$BASE/matrizes/$ID"
  P="$BASE/au/$ID/atuais_Falconidae"
  [ -s "$P/hipoteses.trees" ] || { echo "sem hipoteses: $ID"; continue; }
  grep -v -x -f <(echo morfologia_livre) "$P/ordem.txt" > "$P/ordem_esqueletos.txt"
  N=$(wc -l < "$P/ordem_esqueletos.txt")
  head -n "$N" "$P/hipoteses.trees" > "$P/esqueletos.trees"
  CMD=(iqtree3 -s "$M/$ID.phy" -st MORPH -p "$M/particoes.nex" -te "$BASE/arvores/$ID/esqueleto.treefile"
       -z "$P/esqueletos.trees" -n 0 -zb 10000 -zw -au --seed "$SEMENTE" -T "${THREADS:-1}"
       --prefix "$P/AU_esqueletos" -redo --quiet)
  echo "$(date -Is) ${CMD[*]}" >> "$LOG"
  "${CMD[@]}"
  echo "$ID: $N esqueletos testados"
done
