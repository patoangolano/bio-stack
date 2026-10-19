#!/bin/bash
# Teste do mecanismo da discordancia q x escore bruto (ver research note, secao 3.1).
# So com arvores de genes que tem TODOS os 42 taxa focais (sem dado faltante):
#   (a) colapsadas em aLRT < 0,95 (com politomias);
#   (b) as mesmas arvores sem colapsar (binarias, sem dado faltante).
# Em (b) cada quarteto de especies e resolvido por todas as arvores, entao o q
# (media por quarteto) e o escore bruto (soma) tem de apontar para o mesmo lado.
# Pontua T1 (63K), T2 (Cariamiformes+Falconiformes) e T4 (S+C).
#   bash -l teste_ponderacao.sh <tipo>
set -euo pipefail
O=~/rapinantes/reanalise/b10k/aster
T=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/reanalise/b10k/topologias
tipo=$1
cd "$O"
# arvore com 42 folhas = 41 virgulas (42 campos separados por virgula)
awk -F, 'NF==42' "$tipo.c95.focal.tre" > "$tipo.c95.focal.completas.tre"
awk -F, 'NF==42' "$tipo.alrt.focal.tre" > "$tipo.alrt.focal.completas.tre"
echo "# $tipo: $(wc -l < "$tipo.c95.focal.completas.tre") colapsadas e $(wc -l < "$tipo.alrt.focal.completas.tre") binarias completas"
for conj in c95 alrt; do
  for b in T1_63K T2_Australaves_K+F T4_Afroaves_S+C; do
    astral -C -c "$T/$b.tre" -t 2 --seed 12345 -u 2 -i "$tipo.$conj.focal.completas.tre" -o "completas_${conj}_${tipo}_$b.tre" 2> "completas_${conj}_${tipo}_$b.log"
    printf '%s\t%s\t%s\t%s\n' "$tipo" "$conj" "$b" "$(grep -E '^Score:' "completas_${conj}_${tipo}_$b.log" | tail -1 | awk '{print $2}')"
    cp "completas_${conj}_${tipo}_$b.tre" "$T/"
  done
done
