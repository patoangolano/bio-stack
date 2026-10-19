#!/usr/bin/env bash
# f02e — ortologos NCBI (CDS + proteina + product report) dos primeiros N candidatos
# da permutacao do f02d, filtrados nos taxons RefSeq de taxons.tsv. Lotes de 100 genes.
# Uso: bash f02e_baixar_ortologos.sh <N>
set -euo pipefail
N=${1:-1200}
P=/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica
O=~/rapinantes/locos/ortologos
mkdir -p "$O"; cd "$O"
TAX=$(cat "$P/taxids_refseq.txt")
sed -n "2,$((N+1))p" "$P/candidatos_permutados.tsv" | cut -f2 | split -l 100 -d -a 3 - lote_
for f in lote_[0-9][0-9][0-9]; do
  [ -d "$f.d" ] && { echo "$f ja feito"; continue; }
  for tent in 1 2 3; do
    if datasets download gene gene-id --inputfile "$f" --ortholog "$TAX" --include cds,protein,product-report \
         --filename "$f.zip" --no-progressbar 2>&1 | grep -v "New version"; then :; fi
    if unzip -o -q "$f.zip" -d "$f.tmp" 2>/dev/null; then mv "$f.tmp" "$f.d"; rm -f "$f.zip"; break; fi
    echo "$f tentativa $tent falhou"; rm -rf "$f.tmp" "$f.zip"; sleep 20
  done
  echo "$(date +%T) $f: $(grep -c '>' "$f.d/ncbi_dataset/data/cds.fna" 2>/dev/null) CDS"
done
du -sh "$O"
