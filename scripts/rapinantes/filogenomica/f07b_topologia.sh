#!/usr/bin/env bash
# f07b — buscas ML restritas (uma por hipotese) e teste de topologia (KH, SH, AU)
# na supermatriz de uma execucao do pipeline, sob o esquema de particao que ela
# escolheu (05_concatenado/concatenado.best_model.nex). IQ-TREE 3.1.4 (pixi).
#   bash f07b_topologia.sh <dir da execucao (WSL)> <dir das hipoteses (WSL)> [threads]
# Busca restrita: -g <restricao> -t <inicial por SPR da ML> --fast (a inicial ja esta
# perto do otimo restrito; a ML entra no teste com a topologia da busca completa).
# Hipotese que a ML ja satisfaz nao tem busca: a melhor arvore restrita e a ML.
# Teste: -z (ML, wASTRAL, restritas) -n 0 -zb 10000 -au, semente 12345.
set -euo pipefail
RUN=$1; H=$2; TH=${3:-8}
S="$RUN/05_concatenado/supermatriz.fa"; M="$RUN/05_concatenado/concatenado.best_model.nex"
cd "$H"
iqtree3 --version | head -1 > versao_iqtree.txt
cp ml.nwk candidatos.nwk; echo ML > candidatos.ids
cat astral.nwk >> candidatos.nwk; echo wASTRAL >> candidatos.ids
tail -n +2 hipoteses.tsv | while IFS=$'\t' read -r hid desc clados ja ini; do
  if [ "$ja" = sim ]; then echo "$hid: a ML ja satisfaz; sem busca"; continue; fi
  if [ ! -s "$hid.inicial.nwk" ]; then echo "$hid: sem arvore inicial; busca a partir de parcimonia"; INI=""; else INI="-t $hid.inicial.nwk"; fi
  if [ ! -s "$hid.treefile" ]; then
    echo "$(date +%T) $hid: busca restrita"
    iqtree3 -s "$S" -p "$M" -g "$hid.restricao.nwk" $INI --fast -T "$TH" --seed 12345 --prefix "$hid" -redo --quiet
  fi
  cat "$hid.treefile" >> candidatos.nwk; echo "$hid" >> candidatos.ids
done
echo "$(date +%T) teste de topologia em $(wc -l < candidatos.nwk) arvores"
iqtree3 -s "$S" -p "$M" -t ml.nwk -z candidatos.nwk -n 0 -zb 10000 -au -T "$TH" --seed 12345 --prefix au -redo --quiet
sed -n '/USER TREES/,/^deltaL/p' au.iqtree
