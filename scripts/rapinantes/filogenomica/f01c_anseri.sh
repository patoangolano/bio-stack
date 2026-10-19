#!/usr/bin/env bash
set -euo pipefail
D=~/rapinantes/amostragem
CAMPOS=accession,organism-name,organism-tax-id,assminfo-level,assmstats-contig-n50,assmstats-scaffold-n50,assmstats-total-sequence-len,annotinfo-name,annotinfo-release-date,assminfo-refseq-category,assmstats-gc-percent
datasets summary genome taxon Anseriformes --annotated --assembly-source RefSeq --as-json-lines | dataformat tsv genome --fields "$CAMPOS" > "$D/refseq_Anseriformes.tsv"
cut -f1,2,3,4,5 "$D/refseq_Anseriformes.tsv"
# tax ids e tamanhos das montagens GenBank escolhidas
for a in GCA_034782735.1 GCA_026930845.1 GCA_026212715.1 GCA_028022735.1 GCA_055772235.1 GCA_982185325.1 GCA_054553165.1 GCA_013399355.1 GCA_035590275.1 GCA_037962945.1 GCA_982290615.1; do
  datasets summary genome accession "$a" --as-json-lines | dataformat tsv genome --fields "$CAMPOS" | tail -n +2
done
