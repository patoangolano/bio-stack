#!/usr/bin/env bash
# f02g — substitui o f02e para os lotes que nao baixaram: o servidor do NCBI cortava
# lotes de 100 genes com "stream error ... INTERNAL_ERROR". Mesmos genes, lotes de 20,
# 3 em paralelo, ate 4 tentativas; o gene de um lote que falha 4 vezes vai para
# falhas_lotes.txt (sai do pool por motivo tecnico, nao biologico).
set -uo pipefail
P=/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica
O=~/rapinantes/locos/ortologos; cd "$O"
TAX=$(cat "$P/taxids_refseq.txt"); export TAX
for f in lote_[0-9][0-9][0-9]; do
  [ -d "$f.d" ] && continue
  split -l 20 -d -a 1 "$f" "${f}_p"
done
um() {
  f=$1
  [ -d "$f.d" ] && return 0
  for tent in 1 2 3 4; do
    datasets download gene gene-id --inputfile "$f" --ortholog "$TAX" --include cds,protein,product-report \
      --filename "$f.zip" --no-progressbar > "$f.err" 2>&1
    if unzip -o -q "$f.zip" -d "$f.tmp" 2>/dev/null; then mv "$f.tmp" "$f.d"; rm -f "$f.zip" "$f.err"
      echo "$(date +%T) $f ok ($tent)"; return 0; fi
    rm -rf "$f.tmp" "$f.zip"; sleep 10
  done
  echo "$(date +%T) $f FALHOU"; cat "$f" >> falhas_lotes.txt
}
export -f um
ls lote_[0-9][0-9][0-9]_p[0-9] | grep -v '\.' | xargs -P 3 -I{} bash -c 'um {}'
echo "lotes feitos: $(ls -d lote_*.d | wc -l)"
