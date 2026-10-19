#!/usr/bin/env bash
# Espera a analise livre do LZ07 terminar (f03), interrompe a analise com
# esqueleto que o f03 começaria com UFBoot, e roda a versao sem UFBoot (f06).
set -uo pipefail
BASE=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis
until [ -s "$BASE/arvores/LZ07/livre.ufboot" ] && grep -q "LZ07 livre" "$BASE/log/f03_LZ07.out"; do sleep 20; done
pkill -f "f03_iqtree.sh LZ07" || true
pkill -f "prefix $BASE/arvores/LZ07/esqueleto" || true
sleep 3
rm -f "$BASE/arvores/LZ07/esqueleto."*
THREADS=2 bash /mnt/c/Users/Manec/bio-stack/scripts/rapinantes/fosseis/f06_lz07_esqueleto.sh
