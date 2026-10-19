#!/usr/bin/env bash
# f08a — concordancia so com arvores de genes e wASTRAL (antes da supermatriz):
#   gCF das arvores de genes sobre a topologia do wASTRAL (IQ-TREE 3.1.4, --gcf; o
#   .cf.branch da o ID de cada ramo) e frequencias de quarteto q1/q2/q3 por ramo
#   (wASTRAL --mode 4, nao ponderado, -C -c: so pontua a arvore dada, -u 2).
#   bash f08a_genes_astral.sh <dir da execucao (WSL)> <dir de saida (WSL)>
set -euo pipefail
RUN=$1; O=$2; mkdir -p "$O"; cd "$O"
cat $(ls "$RUN"/04_arvores_genes/*.treefile | sort) > genes.tre
sed -E 's/:nan//g; s/\)[^:,;()]+/)/g; s/:[0-9.eE+-]+//g' "$RUN/07_arvore_especies/especies.tre" > astral_topologia.nwk
iqtree3 -t astral_topologia.nwk --gcf genes.tre --prefix gcf_astral -T 1 --seed 12345 -redo --quiet
wastral -C -c astral_topologia.nwk -i genes.tre --mode 4 -u 2 --seed 12345 -o quartetos_astral.tre 2> quartetos_astral.log
iqtree3 --version | head -1 > versoes.txt; wastral -h 2>&1 | grep Version >> versoes.txt
echo "genes: $(wc -l < genes.tre)"; head -c 600 quartetos_astral.tre; echo
