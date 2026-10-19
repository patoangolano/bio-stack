#!/usr/bin/env bash
# f08d — sCFL (100 quartetos, mesmo modelo particionado e semente da supermatriz) nas arvores restritas do f07b que
# forcam as alternativas NNI dos nos contestados: H2 (Hieraves), H1b (Accipitrimorphae + Coraciimorphae) e H7b
# (Cariama + Falconiformes). O sCF do ramo forcado em cada uma diz qual alternativa da arvore ML e o sDF1 e qual e o
# sDF2 (o .cf.stat nao diz). Leitura: f08e_scfl_alternativas.py.
#   bash f08d_scfl_alternativas.sh <dir da execucao (WSL)> <dir das hipoteses (WSL)> [threads]
set -euo pipefail
RUN=$1; H=$2; TH=${3:-18}
S="$RUN/05_concatenado/supermatriz.fa"; M="$RUN/05_concatenado/concatenado.best_model.nex"
cd "$H"
for h in H2 H1b H7b; do
  iqtree3 -te "$h.treefile" -s "$S" -p "$M" --scfl 100 --prefix "$h.scfl" -T "$TH" --seed 12345 -redo --quiet
done
