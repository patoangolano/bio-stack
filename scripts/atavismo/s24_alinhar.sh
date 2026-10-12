#!/usr/bin/env bash
# Alinhamento do teste s24: os lncRNAs do chrZ do frango sao os restos do
# cluster de esmalte do jacare?
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s24_alinhar.sh
#
# Roda no WSL porque e la que vivem blast e miniprot (camada pixi). Precisa de
# shell de LOGIN: ~/.pixi/bin so entra no PATH em ~/.profile. Sem o -l, blastn
# e miniprot dao command not found.
#
# Insumo: scripts/atavismo/s24_lncrna_relicto.py escreve os FASTA antes.
#
# QUATRO ALINHAMENTOS, nesta ordem, e a ordem importa:
#
#   0. SANIDADE   proteinas do jacare -> regiao do jacare
#                 Se ENAM/AMBN/AMTN nao acertarem o proprio genoma, a janela
#                 esta errada e nada depois vale.
#   1. POSITIVO   LPL do jacare -> intervalo do frango
#                 LPL esta anotado nas duas especies. Recuperar LPL prova que o
#                 metodo tem sensibilidade a 116 Myr. Sem isso, um resultado
#                 negativo sobre o cluster nao se sustenta — seria so falta de
#                 poder de deteccao.
#   2. ALVO       ENAM/AMBN/AMTN do jacare -> intervalo do frango
#   3. RECIPROCO  lncRNAs do frango -> regiao do jacare
#                 Direcao oposta, para nao depender de um unico algoritmo.
set -euo pipefail

PROJ=/mnt/c/Users/Manec/bio-stack
D="$PROJ/runs/atavismo/lncrna"
cd "$D"

echo "=== versoes (vao para Metodos) ==="
tblastn -version | head -1
blastn -version | head -1
miniprot --version 2>&1 | head -1 | sed 's/^/miniprot /'
echo

# Colunas do outfmt 6, nomeadas para o relatorio ficar legivel.
FMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
CAB=$'consulta\tsujeito\t%id\tcompr\tmism\tgaps\tq_ini\tq_fim\ts_ini\ts_fim\tevalue\tbits'

echo "=== bancos ==="
makeblastdb -in frango_chrZ_intervalo.fna -dbtype nucl -out db_frango >/dev/null
makeblastdb -in jacare_chr3_cluster.fna  -dbtype nucl -out db_jacare >/dev/null
echo "  ok"
echo

# -------------------------------------------------------------------------
echo "############ 0. SANIDADE — proteinas do jacare contra o proprio genoma"
echo "############    esperado: os cinco acertam, com identidade alta"
tblastn -query jacare_proteinas.faa -db db_jacare -evalue 1e-5 \
        -seg yes -max_target_seqs 5 -num_threads 4 \
        -outfmt "$FMT" > sanidade.tsv || true
if [ -s sanidade.tsv ]; then
  { echo "$CAB"; sort -k1,1 -k12,12gr sanidade.tsv | awk '!v[$1]++'; } | column -t -s $'\t'
else
  echo "  NENHUM HIT — a janela do jacare esta errada. Pare aqui."
  exit 1
fi
echo

# -------------------------------------------------------------------------
echo "############ 1+2. POSITIVO e ALVO — proteinas do jacare -> chrZ do frango"
echo "############      evalue frouxo (10) de proposito: quero ver hit marginal"
tblastn -query jacare_proteinas.faa -db db_frango -evalue 10 \
        -seg yes -max_target_seqs 20 -num_threads 4 \
        -outfmt "$FMT" > alvo_tblastn.tsv || true
if [ -s alvo_tblastn.tsv ]; then
  { echo "$CAB"; sort -k1,1 -k12,12gr alvo_tblastn.tsv; } | column -t -s $'\t'
else
  echo "  NENHUM HIT para nenhuma das cinco consultas."
fi
echo
echo "--- resumo por consulta (melhor hit de cada) ---"
for sym in LPL ENAM AMBN AMTN DSPP; do
  papel="alvo"
  [ "$sym" = "LPL" ]  && papel="CONTROLE POSITIVO"
  [ "$sym" = "DSPP" ] && papel="CONTROLE NEGATIVO"
  linha=$(awk -v s="${sym}_Amis" '$1==s' alvo_tblastn.tsv 2>/dev/null | sort -k12,12gr | head -1)
  if [ -n "$linha" ]; then
    ev=$(echo "$linha" | cut -f11); bs=$(echo "$linha" | cut -f12)
    pid=$(echo "$linha" | cut -f3); ln=$(echo "$linha" | cut -f4)
    printf '  %-6s %-18s HIT   %%id=%-6s compr=%-5s evalue=%-10s bits=%s\n' \
           "$sym" "$papel" "$pid" "$ln" "$ev" "$bs"
  else
    printf '  %-6s %-18s sem hit\n' "$sym" "$papel"
  fi
done
echo

# -------------------------------------------------------------------------
echo "############ 2b. MINIPROT — alinhamento proteina->genoma com splice"
echo "############      segunda opiniao: o tblastn nao modela intron"
miniprot -t 4 --outs=0.5 --gff frango_chrZ_intervalo.fna jacare_proteinas.faa \
         > alvo_miniprot.gff 2> miniprot.err || true
if grep -qE "^##PAF|\bmRNA\b" alvo_miniprot.gff 2>/dev/null; then
  echo "  consulta | identidade | span no chrZ"
  awk -F'\t' '$3=="mRNA"{
      id="?"; sp="?";
      if (match($9, /Identity=[0-9.]+/))  id=substr($9, RSTART+9, RLENGTH-9);
      if (match($9, /Target=[^;]+/))      sp=substr($9, RSTART+7, RLENGTH-7);
      printf "  %-28s %-10s %s..%s\n", sp, id, $4, $5
  }' alvo_miniprot.gff | sort -u
else
  echo "  miniprot: nenhum alinhamento acima do limiar (--outs=0.5)"
  tail -3 miniprot.err 2>/dev/null | sed 's/^/    /'
fi
echo

# -------------------------------------------------------------------------
echo "############ 3. RECIPROCO — lncRNAs do frango -> regiao do jacare"
echo "############    dc-megablast: desenhado para comparacao entre especies"
blastn -task dc-megablast -query frango_lncrnas.fna -db db_jacare \
       -evalue 10 -dust yes -max_target_seqs 20 -num_threads 4 \
       -outfmt "$FMT" > reciproco_dc.tsv || true
blastn -task blastn -query frango_lncrnas.fna -db db_jacare \
       -evalue 10 -dust yes -word_size 7 -max_target_seqs 20 -num_threads 4 \
       -outfmt "$FMT" > reciproco_sensivel.tsv || true

for arq in reciproco_dc.tsv reciproco_sensivel.tsv; do
  rot=$([ "$arq" = "reciproco_dc.tsv" ] && echo "dc-megablast" || echo "blastn word_size=7")
  echo "--- $rot ---"
  if [ -s "$arq" ]; then
    { echo "$CAB"; sort -k1,1 -k12,12gr "$arq"; } | column -t -s $'\t'
  else
    echo "  nenhum hit"
  fi
done
echo

echo "=== arquivos em $D ==="
ls -la *.tsv *.gff 2>/dev/null | awk '{printf "  %-28s %s bytes\n", $9, $5}'
