#!/usr/bin/env bash
# Verificacao dos hits do s24. Tres testes, e nenhum hit sobrevive sem os tres.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s25_verificar_relicto.sh
#
# O s24 achou sinal em tres lugares. Antes de chamar qualquer coisa de reliquia:
#
#   TESTE A  NUMERO DE COPIAS. Pega a regiao do frango que deu hit e joga contra
#            a janela inteira do jacare. Uma copia = especifico. Muitas = elemento
#            repetitivo, e o hit nao diz nada sobre esmalte. No s24, a consulta
#            ENSGALG00010014761 posicao 3469-3580 bateu em OITO lugares do jacare
#            — assinatura de repeticao, nao de homologia de gene.
#
#   TESTE B  MELHOR HIT RECIPROCO. A regiao do frango, traduzida, tem como melhor
#            acerto a proteina que a encontrou? Se ENAM achou um trecho do chrZ e
#            aquele trecho, de volta, acerta melhor outra proteina, nao e ENAM.
#
#   TESTE C  EXON OU INTRON. O hit nucleotidico do lncRNA cai em exon de AMTN ou
#            em intron? Em intron nao e reliquia de codificante — pode ser
#            elemento nao-codificante conservado ou repeticao compartilhada.
set -euo pipefail

PROJ=/mnt/c/Users/Manec/bio-stack
D="$PROJ/runs/atavismo/lncrna"
cd "$D"

# Janelas, para converter coordenada local em absoluta.
FRANGO_INI=55070000
JACARE_INI=157250000

FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'

echo "=== TESTE A — numero de copias das regioes que deram hit ==="
echo
# Regioes do frango a testar, em coordenada LOCAL da janela (do s24):
#   ENAM  85300-87400   (quatro HSPs agrupados, evalue ate 5e-08)
#   AMBN  98700-99150   (um HSP forte, evalue 1.29e-10)
#   repeticao suspeita: o span do lncRNA 014761 em torno de 3469-3580 ja veio
#   de regiao que bate em oito lugares — entra como controle de repeticao.
cat > regioes.bed <<'BED'
Ggal_chrZ_NRG1_PSD3	85299	87400	ENAM_candidato
Ggal_chrZ_NRG1_PSD3	98699	99150	AMBN_candidato
BED

if ! command -v seqkit >/dev/null; then echo "seqkit ausente"; exit 1; fi
seqkit subseq --bed regioes.bed frango_chrZ_intervalo.fna > regioes_frango.fna 2>/dev/null
echo "regioes extraidas:"
seqkit fx2tab -nl regioes_frango.fna | awk '{printf "  %-50s %s bp\n", $1, $NF}'
echo

echo "--- cada regiao contra a janela do jacare (contagem de locais distintos) ---"
blastn -task dc-megablast -query regioes_frango.fna -db db_jacare \
       -evalue 1e-3 -dust yes -max_target_seqs 50 -num_threads 4 \
       -outfmt "$FMT" > copias_jacare.tsv || true
if [ -s copias_jacare.tsv ]; then
  awk -v ji="$JACARE_INI" -F'\t' '{
      abs_ini = ji + $7 - 1; abs_fim = ji + $8 - 1;
      printf "  %-22s %6.1f%%  %5s bp  jacare:%d-%d  e=%s\n", $1, $3, $4, abs_ini, abs_fim, $9
  }' copias_jacare.tsv
  echo
  echo "  locais distintos por regiao:"
  awk -F'\t' '{print $1}' copias_jacare.tsv | sort | uniq -c | awk '{printf "    %-22s %s hits\n", $2, $1}'
else
  echo "  nenhum hit nucleotidico (esperado: a 116 Myr o nucleotideo satura; o sinal era proteico)"
fi
echo

echo "--- e contra a janela do PROPRIO frango: a regiao e repetitiva em casa? ---"
blastn -task dc-megablast -query regioes_frango.fna -db db_frango \
       -evalue 1e-3 -dust yes -max_target_seqs 50 -num_threads 4 \
       -outfmt "$FMT" > copias_frango.tsv || true
awk -F'\t' '{print $1}' copias_frango.tsv | sort | uniq -c | awk '{printf "    %-22s %s hits na propria janela (1 = unico, >1 = repetido)\n", $2, $1}'
echo

echo "=== TESTE B — melhor hit reciproco das regioes candidatas ==="
echo "    blastx da regiao do frango contra as cinco proteinas do jacare"
echo
makeblastdb -in jacare_proteinas.faa -dbtype prot -out db_prot >/dev/null
blastx -query regioes_frango.fna -db db_prot -evalue 10 -seg yes \
       -max_target_seqs 10 -num_threads 4 \
       -outfmt "$FMT" > reciproco_blastx.tsv || true
if [ -s reciproco_blastx.tsv ]; then
  { printf 'regiao\tproteina\t%%id\tcompr\tq_ini\tq_fim\tp_ini\tp_fim\tevalue\tbits\n'
    sort -k1,1 -k10,10gr reciproco_blastx.tsv; } | column -t -s $'\t'
  echo
  echo "  VEREDITO por regiao (melhor acerto):"
  # O seqkit nomeia a saida pela COORDENADA, nao pelo campo 4 do BED. O
  # casamento e pela coordenada. Bug corrigido em 2026-10-12: a versao anterior
  # filtrava por "ENAM_candidato", nao casava com nada, e imprimia "sem hit"
  # com a tabela cheia de hits logo acima — falha silenciosa do pior tipo,
  # porque o numero errado era zero e zero parece resultado.
  for par in "ENAM:85300-87400" "AMBN:98700-99150"; do
    esperado="${par%%:*}"
    coord="${par#*:}"
    melhor=$(awk -F'\t' -v c="$coord" 'index($1, c)' reciproco_blastx.tsv \
             | sort -t$'\t' -k10,10gr | head -1)
    if [ -n "$melhor" ]; then
      prot=$(printf '%s\n' "$melhor" | cut -f2)
      ev=$(printf  '%s\n' "$melhor" | cut -f9)
      bs=$(printf  '%s\n' "$melhor" | cut -f10)
      n=$(awk -F'\t' -v c="$coord" -v e="${esperado}_Amis" \
          'index($1,c) && $2==e' reciproco_blastx.tsv | wc -l)
      case "$prot" in
        "${esperado}_Amis") v="RECIPROCO CONFIRMA" ;;
        *)                  v="DIVERGE -> melhor acerto e $prot" ;;
      esac
      printf '    %-6s melhor=%-12s e=%-11s bits=%-6s HSPs=%-3s %s\n' \
             "$esperado" "$prot" "$ev" "$bs" "$n" "$v"
    else
      printf '    %-6s sem hit\n' "$esperado"
    fi
  done
else
  echo "  nenhum hit"
fi
echo

echo "=== TESTE C — o hit do lncRNA 014557 cai em exon de AMTN? ==="
echo
# AMTN no jacare: 157.359.292-157.374.521. O hit do s24 foi em 113.281-113.438
# da janela, isto e, 157.363.280-157.363.437 absoluto.
echo "  hit do s24: jacare 157.363.280-157.363.437 (dentro do span de AMTN)"
echo "  AMTN tem 375 aa = ~1.125 bp de CDS espalhados em 15.229 bp de gene."
echo "  Se o hit cair em intron, nao e reliquia de codificante."
echo
echo "  --- a mesma regiao do frango, traduzida, acerta a proteina AMTN? ---"
cat > regiao_amtn.bed <<'BED'
ENSGALG00010014557_Ggal	150	400	lnc557_trecho_do_hit
BED
seqkit subseq --bed regiao_amtn.bed frango_lncrnas.fna > regiao_amtn.fna 2>/dev/null || true
if [ -s regiao_amtn.fna ]; then
  blastx -query regiao_amtn.fna -db db_prot -evalue 100 -seg yes \
         -outfmt "$FMT" -num_threads 4 > amtn_blastx.tsv || true
  if [ -s amtn_blastx.tsv ]; then
    { printf 'regiao\tproteina\t%%id\tcompr\tq_ini\tq_fim\tp_ini\tp_fim\tevalue\tbits\n'
      sort -k10,10gr amtn_blastx.tsv; } | column -t -s $'\t'
  else
    echo "  nenhum hit proteico — o sinal do lncRNA 014557 NAO e codificante."
    echo "  Isso favorece elemento nao-codificante conservado ou repeticao, nao reliquia de CDS."
  fi
else
  echo "  nao consegui extrair o trecho"
fi
echo
echo "=== arquivos ==="
ls -la copias_*.tsv reciproco_blastx.tsv amtn_blastx.tsv 2>/dev/null | awk '{printf "  %-26s %s bytes\n", $9, $5}'
