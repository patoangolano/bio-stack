#!/usr/bin/env bash
# Mede o tempo de uma arvore de gene AA (MFP + UFBoot 1000 --bnni + SH-aLRT 1000 + aBayes, -T 1)
set -euo pipefail
D=~/rapinantes/teste_tempo; mkdir -p "$D"; cd "$D"
for g in PIGT NUBP1; do
  seqkit translate -T 1 --trim ~/rapinantes/locos/refseq_cds/$g.fna 2>/dev/null > $g.faa
  mafft --localpair --maxiterate 1000 --threadit 0 --quiet $g.faa > $g.aln
  /usr/bin/time -f "$g %e s, %M KB" iqtree3 -s $g.aln -m MFP -B 1000 --bnni --alrt 1000 --abayes --seed 12345 -T 1 --prefix $g -redo --quiet || true
  grep -E "Best-fit model|Total wall" $g.iqtree $g.log | head -3
done
