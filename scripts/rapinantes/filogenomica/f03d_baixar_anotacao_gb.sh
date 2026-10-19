#!/usr/bin/env bash
# f03d — proteinas da anotacao NCBI das montagens GenBank de Gypaetus e Sarcoramphus
# (GB_2026_01_30 e GB_2025_12_03): controle do miniprot na mesma montagem.
set -euo pipefail
V=~/rapinantes/validacao; mkdir -p "$V"; cd "$V"
for acc in GCA_028022735.1 GCA_037962945.1; do
  [ -s "$acc.protein.faa" ] && continue
  datasets download genome accession "$acc" --include protein --filename "$acc.zip" --no-progressbar 2>&1 | grep -v "New version" || true
  unzip -o -q "$acc.zip" -d "$acc.x"
  find "$acc.x" -name protein.faa -exec cp {} "$acc.protein.faa" \;
  rm -rf "$acc.x" "$acc.zip"
  echo "$acc: $(grep -c '>' "$acc.protein.faa") proteinas"; head -2 "$acc.protein.faa" | cut -c1-150
done
