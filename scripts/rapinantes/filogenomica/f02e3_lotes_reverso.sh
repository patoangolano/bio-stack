#!/usr/bin/env bash
# f02e3 — baixa, em paralelo ao f02e, os lotes do fim para o comeco (lote_011..lote_008),
# com nomes de arquivo proprios; para no primeiro lote ja feito.
set -uo pipefail
P=/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica
O=~/rapinantes/locos/ortologos; cd "$O"
TAX=$(cat "$P/taxids_refseq.txt")
for f in lote_011 lote_010 lote_009 lote_008; do
  [ -d "$f.d" ] && break
  for tent in 1 2 3 4; do
    datasets download gene gene-id --inputfile "$f" --ortholog "$TAX" --include cds,protein,product-report \
      --filename "$f.r.zip" --no-progressbar 2>&1 | grep -v "New version" | head -1
    if unzip -o -q "$f.r.zip" -d "$f.rtmp" 2>/dev/null; then
      [ -d "$f.d" ] || mv "$f.rtmp" "$f.d"; rm -rf "$f.r.zip" "$f.rtmp"; break
    fi
    rm -rf "$f.rtmp" "$f.r.zip"; sleep 20
  done
  echo "$(date +%T) $f (reverso): $(grep -c '>' "$f.d/ncbi_dataset/data/cds.fna" 2>/dev/null) CDS"
done
