#!/usr/bin/env bash
# Alinhamento comparativo do s26: as proteinas do jacare contra o intervalo
# sintenico NRG1..LPL de cada especie.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s26_alinhar.sh
#
# Mesmos controles do s24, agora por especie:
#   LPL  = controle POSITIVO (anotado em todas; se nao acertar, sem sensibilidade)
#   DSPP = controle NEGATIVO (outro cromossomo no jacare, baixa complexidade)
#   ENAM/AMBN/AMTN = alvo
#
# O jacare entra como controle de sanidade: contra o proprio intervalo, os tres
# alvos tem de acertar perto de 100%.
set -euo pipefail

PROJ=/mnt/c/Users/Manec/bio-stack
D="$PROJ/runs/atavismo/lncrna"
cd "$D"

FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'
PROT=jacare_proteinas.faa

if [ ! -s "$PROT" ]; then
  echo "falta $PROT — rode o s24_lncrna_relicto.py primeiro" >&2
  exit 1
fi

echo "tblastn $(tblastn -version | head -1 | awk '{print $2}')  |  limiar: evalue 10 (frouxo de proposito)"
echo

# rotulo -> descricao, na ordem em que interessa ler
declare -A DESC=(
  [Amis]="Alligator mississippiensis  COM dente  (SANIDADE)"
  [Ggan]="Gavialis gangeticus         COM dente"
  [Acar]="Anolis carolinensis         COM dente"
  [Cpic]="Chrysemys picta             sem dente  (perda independente)"
  [Ggal]="Gallus gallus               sem dente  (Neornithes)"
  [Shab]="Strigops habroptila         sem dente  (Neornithes)"
)
ORDEM=(Amis Ggan Acar Cpic Ggal Shab)

printf '%-6s %-44s %-9s %-9s %-9s %-11s %-9s\n' \
       ROT ESPECIE ENAM AMBN AMTN "LPL(+)" "DSPP(-)"
printf '%s\n' "-------------------------------------------------------------------------------------------------------"

for rot in "${ORDEM[@]}"; do
  fna="intervalo_${rot}.fna"
  [ -s "$fna" ] || { printf '%-6s %-44s  (sem intervalo)\n' "$rot" "${DESC[$rot]}"; continue; }

  db="db_${rot}"
  makeblastdb -in "$fna" -dbtype nucl -out "$db" >/dev/null 2>&1

  tsv="s26_${rot}.tsv"
  tblastn -query "$PROT" -db "$db" -evalue 10 -seg yes \
          -max_target_seqs 50 -num_threads 4 -outfmt "$FMT" > "$tsv" 2>/dev/null || true

  linha=$(printf '%-6s %-44s' "$rot" "${DESC[$rot]}")
  for sym in ENAM AMBN AMTN LPL DSPP; do
    melhor=$(awk -F'\t' -v q="${sym}_Amis" '$1==q' "$tsv" 2>/dev/null \
             | sort -t$'\t' -k10,10gr | head -1)
    if [ -n "$melhor" ]; then
      ev=$(printf '%s\n' "$melhor" | cut -f9)
      linha+=$(printf ' %-9s' "$ev")
    else
      linha+=$(printf ' %-9s' "-")
    fi
  done
  echo "$linha"
done

echo
echo "================================================================================"
echo "DETALHE — melhor HSP por gene e especie, com identidade e comprimento"
echo "================================================================================"
for rot in "${ORDEM[@]}"; do
  tsv="s26_${rot}.tsv"
  [ -s "$tsv" ] || continue
  echo
  echo "--- $rot : ${DESC[$rot]} ---"
  printf '  %-7s %-9s %-8s %-7s %-9s %s\n' GENE EVALUE "%ID" COMPR BITS "HSPs"
  for sym in ENAM AMBN AMTN LPL DSPP; do
    melhor=$(awk -F'\t' -v q="${sym}_Amis" '$1==q' "$tsv" | sort -t$'\t' -k10,10gr | head -1)
    n=$(awk -F'\t' -v q="${sym}_Amis" '$1==q' "$tsv" | wc -l)
    if [ -n "$melhor" ]; then
      printf '  %-7s %-9s %-8s %-7s %-9s %s\n' "$sym" \
        "$(printf '%s\n' "$melhor" | cut -f9)" \
        "$(printf '%s\n' "$melhor" | cut -f3)" \
        "$(printf '%s\n' "$melhor" | cut -f4)" \
        "$(printf '%s\n' "$melhor" | cut -f10)" "$n"
    else
      printf '  %-7s %-9s\n' "$sym" "sem hit"
    fi
  done
done

echo
echo "arquivos: $(ls s26_*.tsv 2>/dev/null | tr '\n' ' ')"
