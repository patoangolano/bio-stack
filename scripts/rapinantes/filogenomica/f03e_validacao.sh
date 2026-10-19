#!/usr/bin/env bash
# f03e — controle do miniprot na mesma montagem: proteina extraida pelo miniprot
# (f03c, registros "ok") de Gypaetus e Sarcoramphus vs. proteoma da anotacao NCBI
# dessas montagens (GB_2026_01_30, GB_2025_12_03), DIAMOND blastp --ultra-sensitive.
# Mede: fracao com melhor hit a >= 98 % de identidade e cobertura da consulta >= 90 %.
set -euo pipefail
V=~/rapinantes/validacao; M=~/rapinantes/miniprot/cds_por_taxon; cd "$V"
diamond version > versao_diamond.txt
for par in GCA_028022735.1:Gypaetus_barbatus GCA_037962945.1:Sarcoramphus_papa; do
  acc=${par%%:*}; tax=${par##*:}
  seqkit translate -T 1 "$M/$tax.fna" > "$tax.mp.faa" 2>/dev/null
  diamond makedb --in "$acc.protein.faa" -d "$acc" --quiet
  diamond blastp -q "$tax.mp.faa" -d "$acc" -o "$tax.vs_anotacao.tsv" --ultra-sensitive -k 1 \
    --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore --threads 4 --quiet
  n=$(grep -c '>' "$tax.mp.faa"); h=$(cut -f1 "$tax.vs_anotacao.tsv" | sort -u | wc -l)
  ok=$(awk -F'\t' '$3>=98 && $7>=90' "$tax.vs_anotacao.tsv" | cut -f1 | sort -u | wc -l)
  idmed=$(awk -F'\t' '{s+=$3} END{printf "%.2f", s/NR}' "$tax.vs_anotacao.tsv")
  echo -e "$tax\t$acc\tregistros=$n\tcom_hit=$h\tident>=98_cob>=90=$ok\tidentidade_media=$idmed"
done | tee validacao_resumo.tsv
