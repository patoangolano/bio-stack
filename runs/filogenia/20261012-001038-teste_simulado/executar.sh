#!/usr/bin/env bash
# bio-stack/filogenia — execucao 20261012-001038-teste_simulado
# gerado por scripts/filogenia/filogenia.py em 2026-10-12T03:10:38+00:00
set -euo pipefail
export NXF_ANSI_LOG=false
mkdir -p ~/nf-work/filogenia/20261012-001038-teste_simulado && cd ~/nf-work/filogenia/20261012-001038-teste_simulado
nextflow run '/mnt/c/Users/Manec/bio-stack/pipelines/filogenia' \
    -profile docker \
    -w ~/nf-work/filogenia/work \
    -resume \
    --input '/mnt/c/Users/Manec/bio-stack/pipelines/filogenia/assets/test/locos' \
    --outdir '/mnt/c/Users/Manec/bio-stack/runs/filogenia/20261012-001038-teste_simulado' \
    --alinhador 'mafft' \
    --limpeza 'clipkit' \
    --outgroup 'Homo sapiens' \
    --nome_taxon 'primeiro' \
    --cpus_concatenado '1' \
    --seed '12345'
