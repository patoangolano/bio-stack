#!/usr/bin/env bash
# Genomas inteiros para as buscas fora do intervalo sintenico (s35/s36):
#   Gallus gallus           GCF_016699485.2  1,05 Gb  ODAPH no genoma inteiro
#   Strigops habroptila     GCF_004027225.2  1,15 Gb  idem, segunda ave (Neoaves)
#   Alligator mississippiensis GCF_030867095.1 2,35 Gb CONTROLE POSITIVO: o mesmo
#                           pipeline tem de achar o ODAPH do jacare (LOC106740521)
#                           como melhor hit do genoma inteiro
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s35_genomas.sh
#
# Ficam em ~/genomas dentro do WSL, nao em /mnt/c: o I/O de ext4 e ordens de
# grandeza mais rapido, e 4,5 GB de genoma nao tem o que fazer no repositorio.
# Idempotente: so baixa o que falta.
set -euo pipefail
G="$HOME/genomas"
mkdir -p "$G"
declare -A ACC=([Ggal]=GCF_016699485.2 [Shab]=GCF_004027225.2 [Amis]=GCF_030867095.1)
for rot in Ggal Shab Amis; do
    fa="$G/${rot}.fna"
    if [ -s "$fa" ]; then
        echo "$rot: ja presente ($(du -h "$fa" | cut -f1))"
        continue
    fi
    zip="$G/${rot}.zip"
    datasets download genome accession "${ACC[$rot]}" --include genome --no-progressbar --filename "$zip"
    unzip -p "$zip" '*/*_genomic.fna' > "$fa"
    rm -f "$zip"
    samtools faidx "$fa"
    echo "$rot: ${ACC[$rot]} $(du -h "$fa" | cut -f1), $(wc -l < "$fa.fai") sequencias"
done
