#!/usr/bin/env bash
# Sonda do ambiente WSL para o braco de filogenomica dos rapinantes.
set -uo pipefail
for t in datasets dataformat miniprot seqkit iqtree3 mafft famsa clipkit trimal wastral macse nextflow docker jq csvtk taxonkit; do
  printf '%-12s %s\n' "$t" "$(command -v "$t" || echo AUSENTE)"
done
datasets --version
miniprot --version
iqtree3 --version | head -1
echo "API_KEY definida: ${NCBI_API_KEY:+sim}"
nproc; free -g | head -2
df -h ~ | tail -1
uptime
ps -eo pid,pcpu,pmem,etime,comm --sort=-pcpu | head -15
ls ~/ | head -40
