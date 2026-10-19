#!/usr/bin/env bash
# bio-stack/filogenia — execucao 20261018-084205-rapinantes_DNA12
# gerado por scripts/filogenia/filogenia.py em 2026-10-18T11:42:05+00:00
set -euo pipefail
export NXF_ANSI_LOG=false
mkdir -p ~/nf-work/filogenia/20261018-084205-rapinantes_DNA12 && cd ~/nf-work/filogenia/20261018-084205-rapinantes_DNA12
nextflow run '/mnt/c/Users/Manec/bio-stack/pipelines/filogenia' \
    -profile docker \
    -w ~/nf-work/filogenia/work \
    -resume \
    --input '/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica/locos_dna12' \
    --outdir '/mnt/c/Users/Manec/bio-stack/runs/filogenia/20261018-084205-rapinantes_DNA12' \
    --alinhador 'nenhum' \
    --limpeza 'nenhuma' \
    --outgroup 'Gallus_gallus' \
    --tipo_seq 'DNA' \
    --nome_taxon 'primeiro' \
    --cpus_concatenado '8' \
    --seed '12345' \
    '--concatenar' 'false' '-c' '/mnt/c/Users/Manec/bio-stack/scripts/rapinantes/filogenomica/limite_2cpus.config'
