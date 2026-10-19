#!/bin/bash
# Quartetos do B10K (Stiller et al. 2024) para os ramos de rapinantes.
#   ETAPA=<etapa> bash -l rodar_aster_b10k.sh <tipo> [<tipo> ...]
#   tipo:  intergenico63K exon_c12 intron uce
#   etapa: tudo (padrao) | focal (= pontua + livre) | pontua | livre | full
# Para cada tipo:
#   1. colapsa ramos com aLRT < 0,95 (como no artigo; o 63K usa o arquivo colapsado oficial);
#   2. poda para os 42 taxa focais (taxa_focais.txt), arvore a arvore;
#   pontua: ASTRAL-IV pontuando (-C) a topologia 63K publicada, podada aos focais (-u 3):
#           os mesmos ramos para todos os tipos de dado, mesmo quando a busca livre diverge;
#   livre:  ASTRAL-IV sem peso, busca livre, nas arvores colapsadas e podadas (-u 3), e
#           wASTRAL hibrido com peso aLRT (-L, --mode 1), busca livre, nas nao colapsadas;
#   full:   ASTRAL-IV pontuando a 63K publicada com os 363 taxa, sem poda (-u 3).
# -u 3 grava q1/q2/q3, f1-f3, pp1-pp3 e EN no rotulo e o freqQuad.csv com as
# quadriparticoes explicitas de cada alternativa (renomeado por execucao).
# Cada etapa pula o que ja tem freqQuad.csv pronto (retomavel depois de queda do WSL).
# Threads: 8 (limite de 10 CPUs combinado para esta analise). Semente fixa.
set -euo pipefail
B=~/rapinantes/reanalise/b10k
P=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/reanalise/b10k
S=/mnt/c/Users/Manec/bio-stack/scripts/rapinantes/reanalise
G=$B/genes
O=$B/aster
T=8
SEED=12345
ETAPA=${ETAPA:-tudo}
mkdir -p "$O"
cd "$B"
# topologia publicada sem suporte nem comprimento (o gotree desta versao nao tem `clear`)
sed -E 's/:[0-9.eE+-]+//g; s/\)[0-9.eE+-]+/)/g' brutos/03_species_trees_63K.tre > "$O/63K_topologia.tre"
python3 -I "$S/podar_arvores.py" "$O/63K_topologia.tre" "$P/taxa_focais.txt" "$O/63K_focal_topologia.tre" > /dev/null

quer() { [ "$ETAPA" = tudo ] || [ "$ETAPA" = "$1" ] || { [ "$ETAPA" = focal ] && [ "$1" != full ]; }; }
feito() { [ -s "$O/$1.freqQuad.csv" ]; }

for tipo in "$@"; do
  echo "== $tipo $ETAPA $(date -Is)"
  if [ "$tipo" = intergenico63K ]; then
    col=$G/intergenico63K_c95.alrt.tre
  else
    col=$G/$tipo.c95.tre
    [ -s "$col" ] || gotree collapse support -s 0.95 -t $T -i "$G/$tipo.alrt.tre" -o "$col"
  fi
  if quer pontua || quer livre; then
    if [ ! -s "$O/$tipo.poda.log" ] || [ "$(wc -l < "$O/$tipo.poda.log")" -lt 2 ]; then
      python3 -I "$S/podar_arvores.py" "$col" "$P/taxa_focais.txt" "$O/$tipo.c95.focal.tre" > "$O/$tipo.poda.log"
      python3 -I "$S/podar_arvores.py" "$G/$tipo.alrt.tre" "$P/taxa_focais.txt" "$O/$tipo.alrt.focal.tre" >> "$O/$tipo.poda.log"
    fi
    cat "$O/$tipo.poda.log"
  fi
  cd "$O"
  if quer pontua && ! feito "$tipo.focal.score63K"; then
    astral -C -c 63K_focal_topologia.tre -t $T --seed $SEED -u 3 --root Gallus_gallus -i "$tipo.c95.focal.tre" -o "$tipo.focal.score63K.tre" 2> "$tipo.focal.score63K.log"
    mv freqQuad.csv "$tipo.focal.score63K.freqQuad.csv"
    echo "pontua ok $(date -Is)"
  fi
  if quer full && ! feito "$tipo.full363.score63K"; then
    astral -C -c 63K_topologia.tre -t $T --seed $SEED -u 3 --root Gallus_gallus -i "$col" -o "$tipo.full363.score63K.tre" 2> "$tipo.full363.score63K.log"
    mv freqQuad.csv "$tipo.full363.score63K.freqQuad.csv"
    echo "full ok $(date -Is)"
  fi
  if quer livre && ! feito "$tipo.focal.astral"; then
    astral -t $T --seed $SEED -u 3 --root Gallus_gallus -i "$tipo.c95.focal.tre" -o "$tipo.focal.astral.tre" 2> "$tipo.focal.astral.log"
    mv freqQuad.csv "$tipo.focal.astral.freqQuad.csv"
    echo "livre astral ok $(date -Is)"
  fi
  if quer livre && ! feito "$tipo.focal.wastral"; then
    wastral -L --mode 1 -t $T --seed $SEED -u 3 --root Gallus_gallus -i "$tipo.alrt.focal.tre" -o "$tipo.focal.wastral.tre" 2> "$tipo.focal.wastral.log"
    mv freqQuad.csv "$tipo.focal.wastral.freqQuad.csv"
    echo "livre wastral ok $(date -Is)"
  fi
  cd "$B"
  echo "fim $tipo $(date -Is)"
done
