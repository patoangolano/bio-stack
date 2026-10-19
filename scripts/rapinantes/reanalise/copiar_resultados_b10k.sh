#!/bin/bash
# Copia os resultados pequenos (arvores, freqQuad, logs, gCF) do WSL para o projeto.
# As arvores de genes (brutas, traduzidas e podadas) ficam no WSL.
set -euo pipefail
B=~/rapinantes/reanalise/b10k
P=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/reanalise/b10k
mkdir -p "$P/aster" "$P/gcf" "$P/brutos_manifesto"
cp "$B"/aster/*.tre "$P/aster/" 2>/dev/null || true
# as arvores de genes podadas (*.focal.tre de entrada) sao grandes: nao copia
find "$P/aster" -name '*.c95.focal.tre' -delete
find "$P/aster" -name '*.alrt.focal.tre' -delete
cp "$B"/aster/*.freqQuad.csv "$B"/aster/*.log "$P/aster/" 2>/dev/null || true
cp "$B"/gcf/*.cf.stat "$B"/gcf/*.cf.branch "$B"/gcf/*.log "$P/gcf/" 2>/dev/null || true
cp "$B"/brutos/manifest_download.tsv "$B"/brutos/sha256.txt "$P/brutos_manifesto/"
cp "$B"/genes/*.mapa.json "$P/brutos_manifesto/"
ls "$P/aster" | wc -l
