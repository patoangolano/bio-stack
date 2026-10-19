#!/usr/bin/env bash
# AU dos atuais no LZ07 sobre as arvores feitas por p04_au.py lz07_spr.
set -euo pipefail
BASE=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis
M="$BASE/matrizes/LZ07"
P="$BASE/au/LZ07/atuais_Falconidae"
CMD=(iqtree3 -s "$M/LZ07.phy" -st MORPH -p "$M/particoes.nex" -te "$BASE/arvores/LZ07/esqueleto.treefile"
     -z "$P/esqueletos.trees" -n 0 -zb 10000 -zw -au --seed 20261017 -T "${THREADS:-2}"
     --prefix "$P/AU_esqueletos" -redo --quiet)
echo "$(date -Is) ${CMD[*]}" >> "$BASE/log/comandos_au.txt"
"${CMD[@]}"
echo "LZ07 AU ok"
