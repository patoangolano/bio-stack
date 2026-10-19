#!/usr/bin/env bash
# LZ07 (Livezey & Zusi 2007; 154 Neornithes x 2325 caracteres variaveis) com
# esqueleto de Stiller et al. 2024, SEM UFBoot: aqui a arvore so serve para
# medir homoplasia (IC/IR) e para o AU dos atuais; os fosseis da matriz nao
# sao rapinantes. As arvores alternativas do AU (Falconidae junto de
# Accipitridae; Cathartidae junto de Ciconiidae) sao feitas depois, por SPR do
# Falconidae/Cathartidae sobre esta arvore, sem nova busca (p04_au.py
# lz07_spr) — conservador: a alternativa nao e reotimizada.
set -euo pipefail
BASE=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis
M="$BASE/matrizes/LZ07"
O="$BASE/arvores/LZ07"
SEMENTE=20261017
T="${THREADS:-2}"
LOG="$BASE/log/comandos_iqtree.txt"
if [ ! -s "$O/esqueleto.treefile" ]; then
  CMD=(iqtree3 -s "$M/LZ07.phy" -st MORPH -p "$M/particoes.nex" -g "$M/esqueleto_stiller.nwk"
       --seed "$SEMENTE" -T "$T" -o Lithornis --prefix "$O/esqueleto" -redo --quiet)
  echo "$(date -Is) ${CMD[*]}" >> "$LOG"
  "${CMD[@]}"
fi
echo "LZ07 esqueleto ok: $(grep -m1 'BEST SCORE' "$O/esqueleto.log" || true)"
