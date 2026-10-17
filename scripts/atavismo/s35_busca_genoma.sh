#!/usr/bin/env bash
# ODAPH e ODAM no GENOMA INTEIRO — o que o s32 nao pode fazer, porque buscou so
# nas janelas sinteticas. Separa "o gene saiu da posicao" de "o gene sumiu".
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s35_busca_genoma.sh
#   uv run python scripts/atavismo/s35_busca_genoma.py
#
# Genomas do s35_genomas.sh em ~/genomas: Amis (controle positivo), Ggal, Shab.
# Consultas:
#   ODAPH  jacare (modelo do s32), Podarcis (ORF do s32), humano, gamba
#   ODAM   humano, Xenopus, Rhinatrema (cecilia: anfibio COM dente)
# CONTROLE POSITIVO: no jacare, o melhor hit de ODAPH do genoma inteiro tem de
# cair no LOC106740521 (NC_081825.1:164.672.404-164.676.187). Se nao cair, o
# metodo nao serve e nenhum zero nas aves vale.
#
# Todo hit fora do jacare com e <= 1e-3 vai a melhor hit reciproco: blastx do
# trecho contra o proteoma inteiro do jacare.
set -euo pipefail
G="$HOME/genomas"
A=/mnt/c/Users/Manec/bio-stack/runs/atavismo
O="$A/genoma"
P="$A/painel/proteoma_Amis/db"
mkdir -p "$O"
cat "$A"/odam_odaph/modelo_odaph_Amis.faa "$A"/odam_odaph/modelo_odaph_Pmur.faa \
    "$A"/odam_odaph/consultas_odaph.faa "$A"/odam_odaph/consultas_odam.faa |
    awk '/^>/{print $1; next} {print}' > "$O/consultas.faa"
echo "consultas: $(grep '^>' "$O/consultas.faa" | tr -d '>' | tr '\n' ' ')"
FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen'

# evalue 10 e decisao por BITSCORE: num genoma de 1-2,3 Gb o e-value de um hit
# real de consulta divergente passa de 1e-2 (medido: so a consulta do proprio
# jacare sobreviveu a 1e-2). E o mesmo motivo de o s26 usar bitscore.
for rot in Amis Ggal Shab; do
    [ -s "$G/db_${rot}.nsq" ] || [ -s "$G/db_${rot}.00.nsq" ] || \
        makeblastdb -in "$G/${rot}.fna" -dbtype nucl -parse_seqids -out "$G/db_${rot}" >/dev/null
    tblastn -query "$O/consultas.faa" -db "$G/db_${rot}" -evalue 10 -seg yes \
            -max_target_seqs 50 -num_threads 16 -outfmt "$FMT" > "$O/s35_${rot}.tsv"
    echo "$rot: $(wc -l < "$O/s35_${rot}.tsv") HSPs (e <= 10)"
done

# reciproco: todo HSP >= 30 bits (acima do piso de ruido medido, 22-30) nas
# aves, e os de ODAM no jacare (e onde se procura a referencia positiva)
: > "$O/s35_reciproco.tsv"
for rot in Amis Ggal Shab; do
    awk -F'\t' -v r="$rot" '$10 >= 30 && (r != "Amis" || $1 ~ /^ODAM/)' "$O/s35_${rot}.tsv" |
    sed 's/ref|\([^|]*\)|/\1/' | sort -t$'\t' -k1,1 -k10,10gr |
    while IFS=$'\t' read -r q s _pid _len _qs _qe ss se ev bits _ql; do
        lo=$(( ss < se ? ss : se )); hi=$(( ss > se ? ss : se ))
        melhor=$(samtools faidx "$G/${rot}.fna" "${s}:${lo}-${hi}" |
                 blastx -db "$P" -max_target_seqs 3 -evalue 10 -outfmt '6 sacc bitscore evalue stitle' 2>/dev/null |
                 head -1 | cut -f2- | cut -c1-120) || melhor=""
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$rot" "$q" "$s" "$lo" "$hi" "$ev" "$bits" "${melhor:-sem hit}" \
            >> "$O/s35_reciproco.tsv"
    done
done
echo "reciproco: $(wc -l < "$O/s35_reciproco.tsv") trechos testados"
