#!/usr/bin/env bash
# Corta o refinamento de UFBoot do LZ07 (estimado em >1 h) e o vigia f07;
# guarda a arvore ML ja encontrada (busca concluida, logL -103014.6) e roda o
# esqueleto com 4 threads.
set -uo pipefail
BASE=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis
pkill -f "f07_troca_lz07.sh" || true
pkill -f "f03_iqtree.sh LZ07" || true
pkill -f "prefix $BASE/arvores/LZ07/livre" || true
sleep 3
cp "$BASE/arvores/LZ07/livre.treefile" "$BASE/arvores/LZ07/livre_ML_sem_ufboot.treefile"
grep -m3 "BEST SCORE\|UPDATE BEST" "$BASE/arvores/LZ07/livre.log" | tail -1
echo "$(date -Is) LZ07 livre: UFBoot interrompido no refinamento (100/1000); arvore ML da busca guardada" >> "$BASE/log/comandos_iqtree.txt"
THREADS=4 bash /mnt/c/Users/Manec/bio-stack/scripts/rapinantes/fosseis/f06_lz07_esqueleto.sh
