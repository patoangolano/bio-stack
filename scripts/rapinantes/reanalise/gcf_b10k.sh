#!/bin/bash
# gCF (Minh, Hahn & Lanfear 2020) das arvores de genes do B10K, podadas aos 42 focais,
# sobre a topologia 63K publicada podada aos focais. Sem alinhamento: so gCF (sCF
# precisaria dos alinhamentos, que nao foram baixados). IQ-TREE 3, 1 thread.
#   bash -l gcf_b10k.sh <tipo> [...]
set -euo pipefail
O=~/rapinantes/reanalise/b10k/aster
G=~/rapinantes/reanalise/b10k/gcf
mkdir -p "$G"
cd "$G"
for tipo in "$@"; do
  iqtree3 -t "$O/63K_focal_topologia.tre" --gcf "$O/$tipo.alrt.focal.tre" --prefix "gcf_$tipo" -T 1 --seed 12345 > "gcf_$tipo.stdout" 2>&1
  echo "gcf $tipo ok $(date -Is)"
done
