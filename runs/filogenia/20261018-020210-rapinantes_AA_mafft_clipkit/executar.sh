#!/usr/bin/env bash
# bio-stack/filogenia — execucao 20261018-020210-rapinantes_AA_mafft_clipkit
# gerado por scripts/filogenia/filogenia.py em 2026-10-18T05:02:10+00:00
set -euo pipefail
export NXF_ANSI_LOG=false
mkdir -p ~/nf-work/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit && cd ~/nf-work/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit
nextflow run '/mnt/c/Users/Manec/bio-stack/pipelines/filogenia' \
    -profile docker \
    -w ~/nf-work/filogenia/work \
    -resume \
    --input '/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica/locos_aa' \
    --outdir '/mnt/c/Users/Manec/bio-stack/runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit' \
    --alinhador 'mafft' \
    --limpeza 'clipkit' \
    --outgroup 'Gallus_gallus' \
    --tipo_seq 'AA' \
    --nome_taxon 'primeiro' \
    --cpus_concatenado '8' \
    --seed '12345' \
    '-c' '/mnt/c/Users/Manec/bio-stack/scripts/rapinantes/filogenomica/limite_10cpus.config'
