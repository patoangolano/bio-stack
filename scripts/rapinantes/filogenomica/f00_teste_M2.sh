#!/usr/bin/env bash
# compara miniprot -M 1 (padrao) e -M 2 (metade da amostragem de modimizadores, menos memoria) no Sagittarius
set -euo pipefail
O=~/rapinantes/miniprot_teste; mkdir -p "$O"; cd "$O"
acc=GCA_034782735.1; tmp=~/rapinantes/raw/$acc/genoma.fna
gzip -dc ~/rapinantes/raw/$acc/$acc.fna.gz > "$tmp"
/usr/bin/time -f "M2: %e s, %M KB" miniprot -M 2 -t 8 --gff -P Sagittarius_serpentarius "$tmp" ~/rapinantes/locos/consultas/Aquila_chrysaetos.faa > M2.gff 2> M2.log
tail -1 M2.log
rm -f "$tmp"
for f in ~/rapinantes/miniprot/Sagittarius_serpentarius.gff M2.gff; do
  grep -P "\tmRNA\t" $f | grep "Rank=1;" | awk -F'\t' '{split($9,a,";"); t=""; for(i in a){if(a[i]~/^Target=/)t=a[i]}; print t"\t"$1"\t"$4"\t"$5}' | sed 's/Target=//' | sort > $(basename $f).r1.tsv
done
echo "Rank1 M1: $(wc -l < Sagittarius_serpentarius.gff.r1.tsv); M2: $(wc -l < M2.gff.r1.tsv); identicos (loco, contig, inicio, fim): $(comm -12 Sagittarius_serpentarius.gff.r1.tsv M2.gff.r1.tsv | wc -l)"
