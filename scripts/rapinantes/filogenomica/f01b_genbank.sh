#!/usr/bin/env bash
# f01b — montagens GenBank de Accipitridae e Strigidae (para subfamilias sem RefSeq)
set -euo pipefail
D=~/rapinantes/amostragem
CAMPOS=accession,organism-name,organism-tax-id,assminfo-level,assmstats-contig-n50,assmstats-scaffold-n50,assmstats-total-sequence-len,annotinfo-name,assminfo-release-date
for c in Accipitridae Strigidae; do
  datasets summary genome taxon "$c" --as-json-lines 2>/dev/null | dataformat tsv genome --fields "$CAMPOS" > "$D/genbank_$c.tsv" || true
  echo "$c: $(($(wc -l < "$D/genbank_$c.tsv") - 1))"
done
