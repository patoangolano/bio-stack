set -euo pipefail
O=~/rapinantes/miniprot_teste; cd "$O"
acc=GCA_034782735.1; tmp=~/rapinantes/raw/$acc/genoma.fna
gzip -dc ~/rapinantes/raw/$acc/$acc.fna.gz > "$tmp"
grep -v "^##PAF" Sagittarius_serpentarius.gff > s.gff
gffread -g "$tmp" -x Sagittarius_serpentarius.cds.fna s.gff 2> gffread.log || cat gffread.log
rm -f "$tmp" "$tmp.fai"
grep -c '>' Sagittarius_serpentarius.cds.fna; head -2 Sagittarius_serpentarius.cds.fna | cut -c1-100
mkdir -p ~/rapinantes/miniprot && cp Sagittarius_serpentarius.gff Sagittarius_serpentarius.cds.fna ~/rapinantes/miniprot/
