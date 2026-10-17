#!/usr/bin/env bash
# Alinhamento do s27: para cada gene do painel, TODAS as proteinas do intervalo
# do jacare (alvo, as duas ancoras, os competidores e o LPL negativo) contra o
# intervalo sintenico de cada especie.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s27_alinhar.sh
#
# Mesmo tblastn e mesmos parametros do s26 (evalue 10, SEG ligado), para que os
# bitscores dos dois scripts sejam comparaveis. A decisao de quem fica com cada
# trecho do genoma — alvo ou paralogo — nao e feita aqui: e do s27_indice.py,
# sobre a tabela completa.
set -euo pipefail

D=/mnt/c/Users/Manec/bio-stack/runs/atavismo/painel
cd "$D"

FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'
echo "tblastn $(tblastn -version | head -1 | awk '{print $2}')  |  evalue 10, SEG ligado"

for faa in painel_*_proteinas.faa; do
    alvo=${faa#painel_}
    alvo=${alvo%_proteinas.faa}
    nq=$(grep -c '^>' "$faa")
    echo
    echo "== $alvo  ($nq consultas) =="
    for fna in intervalo_"$alvo"_*.fna; do
        [ -s "$fna" ] || continue
        rot=${fna#intervalo_"$alvo"_}
        rot=${rot%.fna}
        db="db_${alvo}_${rot}"
        makeblastdb -in "$fna" -dbtype nucl -out "$db" >/dev/null
        tsv="s27_${alvo}_${rot}.tsv"
        tblastn -query "$faa" -db "$db" -evalue 10 -seg yes \
                -max_target_seqs 50 -num_threads 8 -outfmt "$FMT" > "$tsv"
        printf '   %-5s %6s HSPs  %s\n' "$rot" "$(wc -l < "$tsv")" "$tsv"
    done
done
