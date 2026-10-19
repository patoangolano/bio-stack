#!/usr/bin/env bash
# f02b — teste do download de ortologos (3 genes do galo) filtrado nos 46 taxons RefSeq
set -euo pipefail
D=~/rapinantes/locos/teste
rm -rf "$D"; mkdir -p "$D"; cd "$D"
TAX=$(cat /mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica/taxids_refseq.txt)
time datasets download gene gene-id 770363 395456 396386 --ortholog "$TAX" --include cds,protein,product-report --filename t.zip --no-progressbar
unzip -o -q t.zip -d t
find t -type f | xargs ls -la
grep -c '>' t/ncbi_dataset/data/cds.fna t/ncbi_dataset/data/protein.faa
head -3 t/ncbi_dataset/data/cds.fna
head -c 1500 t/ncbi_dataset/data/product_report.jsonl; echo
