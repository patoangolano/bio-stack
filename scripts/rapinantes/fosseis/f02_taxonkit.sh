#!/usr/bin/env bash
# Ordem e familia (NCBI Taxonomy, taxdump local em ~/.taxonkit) para cada nome
# de runs/rapinantes/fosseis/taxonomia/nomes.txt. Rodar com bash -l (pixi no PATH).
set -euo pipefail
D=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/fosseis/taxonomia
taxonkit version
cat ~/.taxonkit/BAIXADO_EM
taxonkit name2taxid "$D/nomes.txt" \
  | taxonkit reformat -I 2 -f "{o}\t{f}" \
  > "$D/linhagem.tsv"
wc -l "$D/linhagem.tsv"
