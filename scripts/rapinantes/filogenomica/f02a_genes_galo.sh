#!/usr/bin/env bash
# f02a — relatorio de genes do Gallus gallus (NCBI Gene), para o sorteio de locos.
set -euo pipefail
D=~/rapinantes/locos
mkdir -p "$D"
datasets summary gene taxon 9031 --as-json-lines > "$D/galo_genes.jsonl"
wc -l "$D/galo_genes.jsonl"; ls -la "$D/galo_genes.jsonl"
head -c 3000 "$D/galo_genes.jsonl"; echo
