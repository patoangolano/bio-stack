#!/usr/bin/env bash
set -euo pipefail
cd /mnt/c/Users/user/bio-stack
python3 scripts/energia_microbiana/atualizar_genomas.py \
  --candidates scripts/energia_microbiana/candidatos.tsv \
  --out runs/energia_microbiana_2026-11-17/genomas
