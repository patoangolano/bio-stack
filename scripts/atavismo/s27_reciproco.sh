#!/usr/bin/env bash
# Melhor hit reciproco dos HSPs que sobreviveram a competicao do s27_indice.py,
# contra o PROTEOMA INTEIRO do jacare (GCF_030867095.1, rAllMis1).
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s27_reciproco.sh
#
# Por que o proteoma inteiro: a competicao do s27_indice so enxerga os genes do
# intervalo do jacare. Um paralogo cujo ortologo do jacare mora em outro lugar
# (ou nao esta anotado ali) nao concorre, e o trecho dele pode ficar com o alvo
# por falta de rival. Aqui concorrem as ~40 mil proteinas do jacare.
#
# Veredito por trecho: RECIPROCO se o melhor hit do blastx e o mesmo produto do
# alvo (isoformas contam como o mesmo); senao, nomeia quem ganhou.
# So entram trechos com >= 40 bits: o LPL, controle negativo, da 22-30 bits em
# todos os intervalos, e abaixo disso nao ha o que testar.
set -euo pipefail

D=/mnt/c/Users/Manec/bio-stack/runs/atavismo/painel
P="$D/proteoma_Amis"
cd "$D"

# --- proteoma do jacare -------------------------------------------------------------
if [ ! -s "$P/protein.faa" ]; then
    mkdir -p "$P"
    datasets download genome accession GCF_030867095.1 --include protein \
        --no-progressbar --filename "$P/proteoma.zip" >/dev/null
    unzip -o -q -j "$P/proteoma.zip" '*/protein.faa' -d "$P"
    rm -f "$P/proteoma.zip"
fi
[ -s "$P/db.pin" ] || makeblastdb -in "$P/protein.faa" -dbtype prot -out "$P/db" >/dev/null
echo "proteoma do jacare: $(grep -c '^>' "$P/protein.faa") proteinas (GCF_030867095.1)"

# Produto do alvo = titulo NCBI da proteina consulta, gravado pelo s27_indice.py.
# Nao pelo acesso: a tabela de features cita versoes (XP_014463267.1) que o
# proteoma atual ja substituiu, e um grep vazio aqui derrubava o laco em silencio.
produto() { awk -F'\t' -v a="$1" '$1 == a {gsub(/\r/, "", $2); print $2}' s27_produtos.tsv; }

rm -f s27_reciproco.tmp
echo
printf '%-6s %-5s %-19s %6s   %-10s %s\n' ALVO ESP TRECHO BITS VEREDITO "MELHOR HIT NO PROTEOMA (bits)"
printf '%s\n' "--------------------------------------------------------------------------------------------------"
tail -n +2 s27_sobreviventes.tsv | awk -F'\t' '$6 >= 40' | sort -t$'\t' -k1,1 -k2,2 -k6,6gr |
while IFS=$'\t' read -r alvo esp arq ini fim bits _ev _lpl; do
    id=$(head -1 "$arq" | cut -c2- | awk '{print $1}')
    [ -s "$arq.fai" ] || samtools faidx "$arq"
    esperado=$(produto "$alvo")
    melhor=$(samtools faidx "$arq" "${id}:${ini}-${fim}" |
             blastx -db "$P/db" -max_target_seqs 5 -evalue 10 \
                    -outfmt '6 sacc bitscore stitle' 2>/dev/null | head -1) || melhor=""
    if [ -z "$melhor" ]; then
        veredito="SEM HIT"; nome="-"
    else
        # stitle vem como "XP_... produto [organismo]": tira o acesso e o organismo
        nome=$(printf '%s' "$melhor" | cut -f3 | cut -d' ' -f2- | sed -E 's/ isoform X[0-9]+//; s/ \[.*//')
        b=$(printf '%s' "$melhor" | cut -f2)
        if [ "$nome" = "$esperado" ]; then veredito="RECIPROCO"; else veredito="OUTRO"; fi
        nome="$nome ($b)"
    fi
    printf '%-6s %-5s %-19s %6s   %-10s %s\n' "$alvo" "$esp" "${ini}-${fim}" "$bits" "$veredito" "$nome"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$alvo" "$esp" "$ini" "$fim" "$bits" "$veredito" "$nome" >> s27_reciproco.tmp
done

# --- resumo por gene x especie: so o que passou no reciproco ---------------------------
mv s27_reciproco.tmp s27_reciproco.tsv
echo
echo "RESUMO — trechos confirmados (melhor hit reciproco = o proprio alvo)"
printf '%-6s %-5s %10s %14s %12s\n' ALVO ESP CONFIRMADOS "BITS(soma)" REJEITADOS
awk -F'\t' '{k=$1" "$2; if ($6=="RECIPROCO") {n[k]++; s[k]+=$5} else r[k]++; v[k]=1}
     END {for (k in v) {split(k,a," "); printf "%-6s %-5s %10d %14.1f %12d\n", a[1], a[2], n[k], s[k], r[k]}}' \
    s27_reciproco.tsv | sort
