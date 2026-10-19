#!/usr/bin/env bash
# f01 — amostragem de taxons: genomas RefSeq anotados por clado (NCBI datasets 18.38)
# e, para clados essenciais sem RefSeq, as montagens GenBank disponiveis.
# Saida: ~/rapinantes/amostragem/{refseq,genbank}_<clado>.tsv
set -euo pipefail
D=~/rapinantes/amostragem
mkdir -p "$D"
CAMPOS=accession,organism-name,organism-tax-id,assminfo-level,assmstats-contig-n50,assmstats-scaffold-n50,assmstats-total-sequence-len,annotinfo-name,annotinfo-release-date,assminfo-refseq-category,assmstats-gc-percent
CLADOS="Cathartiformes Accipitriformes Falconiformes Strigiformes Cariamiformes Psittaciformes Passeriformes Coliiformes Leptosomiformes Trogoniformes Bucerotiformes Coraciiformes Piciformes Caprimulgiformes Apodiformes Ciconiiformes Charadriiformes Galliformes Columbiformes Pelecaniformes Opisthocomiformes Gruiformes Otidiformes Musophagiformes Cuculiformes Phaethontiformes Eurypygiformes Gaviiformes Sphenisciformes Procellariiformes Suliformes"
for c in $CLADOS; do
  datasets summary genome taxon "$c" --annotated --assembly-source RefSeq --as-json-lines 2>/dev/null \
    | dataformat tsv genome --fields "$CAMPOS" > "$D/refseq_$c.tsv" || true
  n=$(($(wc -l < "$D/refseq_$c.tsv") - 1))
  echo "$c RefSeq anotados: $n"
done
# GenBank (sem filtro de anotacao) para os clados que provavelmente faltam no RefSeq
for c in Sagittariidae Pandionidae Cathartidae Cariamiformes Leptosomiformes Coliiformes Trogoniformes Falconidae Tytonidae Caprimulgidae Ciconiiformes Strigopidae Acanthisittidae; do
  datasets summary genome taxon "$c" --as-json-lines 2>/dev/null \
    | dataformat tsv genome --fields "$CAMPOS" > "$D/genbank_$c.tsv" || true
  n=$(($(wc -l < "$D/genbank_$c.tsv") - 1))
  echo "$c todas as montagens: $n"
done
