#!/usr/bin/env bash
# f02e2 — refaz lotes do f02e que falharam (erro HTTP/2 do servidor): divide em
# sublotes de 10 genes e, se o sublote falhar 3 vezes, gene a gene. Gene que falha
# sozinho fica registrado em falhas.txt (sai do pool; nao e criterio biologico).
set -uo pipefail
P=/mnt/c/Users/Manec/bio-stack/data/proc/rapinantes/filogenomica
O=~/rapinantes/locos/ortologos
cd "$O"
TAX=$(cat "$P/taxids_refseq.txt")
baixa() {  # $1 = arquivo de ids, $2 = destino .d
  for tent in 1 2 3; do
    datasets download gene gene-id --inputfile "$1" --ortholog "$TAX" --include cds,protein,product-report \
      --filename "$1.zip" --no-progressbar 2>&1 | grep -v "New version" | head -2
    if unzip -o -q "$1.zip" -d "$2.tmp" 2>/dev/null; then mv "$2.tmp" "$2"; rm -f "$1.zip"; return 0; fi
    rm -rf "$2.tmp" "$1.zip"; sleep 15
  done
  return 1
}
for f in lote_[0-9][0-9][0-9]; do
  [ -d "$f.d" ] && continue
  echo "$(date +%T) refazendo $f"
  split -l 10 -d -a 2 "$f" "${f}_s"
  for s in ${f}_s[0-9][0-9]; do
    [ -d "$s.d" ] && continue
    if ! baixa "$s" "$s.d"; then
      while read -r gid; do
        echo "$gid" > "${s}_g$gid"
        baixa "${s}_g$gid" "${s}_g$gid.d" || echo "$gid" >> falhas.txt
      done < "$s"
    fi
  done
done
ls -d lote_*.d | wc -l; [ -f falhas.txt ] && wc -l falhas.txt
