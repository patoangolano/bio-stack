#!/usr/bin/env bash
# teste do miniprot (Sagittarius) para depurar o f03c antes da rodada completa
set -euo pipefail
O=~/rapinantes/miniprot_teste; mkdir -p "$O"; cd "$O"
acc=GCA_034782735.1; tmp=~/rapinantes/raw/$acc/genoma.fna
gzip -dc ~/rapinantes/raw/$acc/$acc.fna.gz > "$tmp"
/usr/bin/time -f "%e s, %M KB" miniprot -t 8 --gff -P Sagittarius_serpentarius "$tmp" ~/rapinantes/locos/consultas/Aquila_chrysaetos.faa > Sagittarius_serpentarius.gff 2> mp.log
tail -1 mp.log
gffread -g "$tmp" -x Sagittarius_serpentarius.cds.fna Sagittarius_serpentarius.gff 2> gffread.log
rm -f "$tmp" "$tmp.fai"
grep -c $'\tmRNA\t' Sagittarius_serpentarius.gff; grep -c '>' Sagittarius_serpentarius.cds.fna
grep $'\tmRNA\t' Sagittarius_serpentarius.gff | head -3
grep -v "^#" Sagittarius_serpentarius.gff | cut -f3 | sort | uniq -c
