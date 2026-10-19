#!/usr/bin/env bash
# f02c — quais tax IDs tem dados no NCBI Orthologs (gene BSG do galo, 770363, como sonda)
set -uo pipefail
T=/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica/taxons.tsv
tail -n +2 "$T" | while IFS=$'\t' read -r taxon esp acc fonte tax resto; do
  [ "$fonte" = RefSeq ] || continue
  out=$(datasets summary gene gene-id 770363 --ortholog "$tax" --as-json-lines 2>&1 | grep -v "New version")
  if echo "$out" | grep -q "no gene orthologs"; then st=SEM_ORTOLOGOS
  elif echo "$out" | grep -q '"gene_id"'; then st="ok:$(echo "$out" | grep -c '"gene_id"')"
  else st="vazio"; fi
  echo -e "$taxon\t$tax\t$st"
done
