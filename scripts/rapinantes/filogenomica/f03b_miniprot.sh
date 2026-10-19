#!/usr/bin/env bash
# f03b — miniprot 0.18 nas montagens GenBank (taxons "miniprot" de taxons.tsv): proteina
# do parente RefSeq (consulta_miniprot) -> genoma, GFF3 e CDS extraido com gffread.
# Uma montagem por vez: descompacta, alinha, extrai e apaga o .fna descompactado.
# Uso: bash f03b_miniprot.sh [threads]   (padrao 8; o projeto limita a 10 CPUs)
set -euo pipefail
TH=${1:-8}
T=/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica/taxons.tsv
R=~/rapinantes/raw; C=~/rapinantes/locos/consultas; O=~/rapinantes/miniprot
mkdir -p "$O"
miniprot --version > "$O/versao_miniprot.txt"; gffread --version > "$O/versao_gffread.txt"
tail -n +2 "$T" | awk -F'\t' '$4=="miniprot"{print $3"\t"$1"\t"$12}' | while IFS=$'\t' read -r acc taxon consulta; do
  g="$R/$acc/$acc.fna.gz"
  if [ -s "$O/$taxon.cds.fna" ]; then echo "$taxon ja feito"; continue; fi
  if [ ! -s "$g" ]; then echo "$taxon: genoma $g ausente"; continue; fi
  echo "$(date +%T) $taxon ($acc) <- $consulta"
  tmp=~/rapinantes/tmp/$acc.fna; mkdir -p ~/rapinantes/tmp   # descompactado fora de raw/
  gzip -dc "$g" > "$tmp"
  # -M 2: metade da amostragem de modimizadores no indice. Medido no Sagittarius: pico de
  # 10,2 GB (-M 1) para 6,5 GB, e 684 de 686 alinhamentos Rank=1 no mesmo lugar. Com
  # -M 1 o WSL (15 GB, dividido com outra analise) reiniciou durante o Pandion.
  /usr/bin/time -v miniprot -M 2 -t "$TH" --gff -P "$taxon" "$tmp" "$C/$consulta.faa" > "$O/$taxon.gff" 2> "$O/$taxon.miniprot.log"
  grep -E "Maximum resident|Elapsed" "$O/$taxon.miniprot.log" || true
  # gffread nao aceita as linhas ##PAF que o miniprot intercala no GFF3
  grep -v "^##PAF" "$O/$taxon.gff" > "$O/$taxon.semPAF.gff"
  gffread -g "$tmp" -x "$O/$taxon.cds.fna" "$O/$taxon.semPAF.gff" 2> "$O/$taxon.gffread.log"
  rm -f "$tmp" "$tmp.fai"
  echo "  mRNA: $(grep -c $'\tmRNA\t' "$O/$taxon.gff"); CDS extraidos: $(grep -c '>' "$O/$taxon.cds.fna")"
done
