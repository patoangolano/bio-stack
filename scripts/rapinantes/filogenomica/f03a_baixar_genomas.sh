#!/usr/bin/env bash
# f03a — baixa as montagens GenBank (so a sequencia) dos taxons marcados "miniprot"
# em taxons.tsv para ~/rapinantes/raw/<acessao>/. Teto do projeto: ~10-11 genomas, <= 15 GB.
set -euo pipefail
T=/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica/taxons.tsv
R=~/rapinantes/raw
mkdir -p "$R"
tail -n +2 "$T" | awk -F'\t' '$4=="miniprot"{print $3"\t"$1}' | while IFS=$'\t' read -r acc taxon; do
  d="$R/$acc"; mkdir -p "$d"
  if ls "$d"/*.fna.gz >/dev/null 2>&1; then echo "$taxon $acc ja baixado"; continue; fi
  echo "$(date +%T) $taxon $acc"
  datasets download genome accession "$acc" --include genome --filename "$d/pacote.zip" --no-progressbar 2>&1 | grep -v "New version" || true
  unzip -o -q "$d/pacote.zip" -d "$d/x"
  f=$(find "$d/x" -name '*.fna' | head -1)
  gzip -c "$f" > "$d/$acc.fna.gz"
  rm -rf "$d/x" "$d/pacote.zip"
  ls -la "$d/$acc.fna.gz"
done
du -sh "$R"
