#!/usr/bin/env bash
# Baixa o MAF do multiz100way do hg38 (chr4, 4,6 GB) e recorta os exons
# codificantes de ENAM/AMBN/AMTN (cds_hg38.bed, do s36_multiz.py bed).
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s36_multiz.sh
#
# O MAF fica em ~/genomas/multiz (ext4); so o recorte (KB) vai para o projeto.
set -euo pipefail
M="$HOME/genomas/multiz"
O=/mnt/c/Users/Manec/bio-stack/runs/atavismo/multiz
URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/chr4.maf.gz
mkdir -p "$M"
if [ ! -s "$M/chr4.maf.gz" ] || [ "$(stat -c %s "$M/chr4.maf.gz")" -ne 4627971825 ]; then
    curl -fsSL --retry 5 -C - -o "$M/chr4.maf.gz" "$URL"
fi
echo "chr4.maf.gz: $(stat -c %s "$M/chr4.maf.gz") bytes"
# as ferramentas kent abrem .gz direto (lineFile descomprime em fluxo); o MAF
# descomprimido teria dezenas de GB e nao precisa existir em disco
mafsInRegion "$O/cds_hg38.bed" "$O/cds_hg38.maf" "$M/chr4.maf.gz"
echo "recorte: $(grep -c '^a' "$O/cds_hg38.maf") blocos em $O/cds_hg38.maf"
