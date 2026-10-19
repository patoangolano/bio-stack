#!/usr/bin/env bash
# pos_concat_020210.sh — the 020210 steps the SIGHUP of 2026-10-18 21:09 cut off, replicated outside Nextflow with the
# pipeline's own containers, arguments and publish patterns: IQTREE_CONCAT publish, IQTREE_CONCORDANCIA, RESUMO.
# Nextflow cannot -resume them: the killed IQTREE_CONCAT never reached its cache entry. The supermatrix itself is the
# fixed-model rerun in $RD (roda_fixo.sh); inputs still come from the original work dir $TD.
set -euo pipefail
TD=/home/manec/nf-work/filogenia/work/81/e3fdfb681f0da353c2d65ff89f587d
RD=${RD:-/home/manec/nf-work/filogenia/manual/020210/fixo/run}
R=/mnt/c/Users/Manec/bio-stack/runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit
BIN=/mnt/c/Users/Manec/bio-stack/pipelines/filogenia/bin
M=/home/manec/nf-work/filogenia/manual/020210
IMG_IQ=quay.io/biocontainers/iqtree:3.1.4--h8471819_0
IMG_PY=quay.io/biocontainers/python:3.12.12@sha256:4cc84261c4b7a77a23cfbaa5b5316416a507b42ac8cde330b002631666d99ca2
SEED=12345
SCFL=$(grep -oE '"scfl": *[0-9]+' "$R/pipeline_info/parametros.json" | grep -oE '[0-9]+$')
CONCAT_DIR=$(dirname "$(readlink -f "$TD/supermatriz.fa")")
docker_em() { local dir=$1 img=$2; shift 2
  docker run --rm -u "$(id -u):$(id -g)" -v /home/manec/nf-work:/home/manec/nf-work -w "$dir" "$img" /bin/bash -euc "$*"; }

echo "[pos] 1/3 publica IQTREE_CONCAT ($(date '+%T'))"
mkdir -p "$R/05_concatenado"
for f in "$RD"/concatenado.{treefile,contree,iqtree,log,best_scheme.nex,best_model.nex}; do
  if [ -f "$f" ]; then cp -p "$f" "$R/05_concatenado/"; fi
done

echo "[pos] 2/3 IQTREE_CONCORDANCIA, sCFL com $SCFL quartetos ($(date '+%T'))"
C=$M/concordancia
rm -rf "$C"; mkdir -p "$C/genes"
cp -p "$RD/concatenado.treefile" "$RD/concatenado.best_model.nex" "$C/"
ln -s "$CONCAT_DIR/supermatriz.fa" "$C/supermatriz.fa"
cp -p "$R"/04_arvores_genes/*.treefile "$C/genes/"
# Same commands as modules/local/iqtree_concordancia/main.nf. The module runs with cpus = 1; sCFL here gets 8 threads
# (approved 2026-10-19 to shorten the run), so its values need not be bit-identical to a -T 1 run. gCF only counts trees.
docker_em "$C" "$IMG_IQ" "
cat genes/*.treefile > arvores_genes.tre
iqtree -t concatenado.treefile --gcf arvores_genes.tre --prefix concordancia.gcf -T 1 --seed $SEED
iqtree -te concatenado.treefile -s supermatriz.fa -p concatenado.best_model.nex --scfl $SCFL --prefix concordancia.scfl -T 8 --seed $SEED
iqtree -version 2>&1 | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -n1 > versao_iqtree.txt"
mkdir -p "$R/06_concordancia"
cp -p "$C"/concordancia.{gcf,scfl}.cf.tree "$C"/concordancia.{gcf,scfl}.cf.stat "$C"/concordancia.{gcf,scfl}.log "$R/06_concordancia/"

echo "[pos] 3/3 RESUMO ($(date '+%T'))"
S=$M/resumo
rm -rf "$S"; mkdir -p "$S"/{preparos,alinhados,limpos,genes,concat,cf,especies}
cp -p "$R/pipeline_info/parametros.json" "$BIN/filogenia_util.py" "$S/"
cp -p "$R"/01_entrada/*.entrada.json "$S/preparos/"
cp -p "$R"/02_alinhamento/*.fas "$S/alinhados/"
cp -p "$R"/03_limpeza/*.limpo.fa "$S/limpos/"
cp -p "$R"/04_arvores_genes/*.iqtree "$R"/04_arvores_genes/*.treefile "$S/genes/"
cp -p "$RD"/concatenado.{iqtree,treefile,log} "$RD"/concatenado*.nex "$CONCAT_DIR/particoes.nex" "$CONCAT_DIR/ocupacao.tsv" "$S/concat/"
cp -p "$C"/concordancia.{gcf,scfl}.cf.stat "$S/cf/"
cp -p "$R/07_arvore_especies/especies.tre" "$R/07_arvore_especies/especies.log" "$S/especies/"

# versoes.yml as the pipeline's collectFile would have written it (it never did: the SIGHUP came first).
# Each value is what that process captured at run time in its .command.env; see versoes_020210.sh.
v_cf=$(cat "$C/versao_iqtree.txt")
# roda_fixo.sh records the version of the image that ran the supermatrix; same image as the concordance if missing.
v_concat=$(cat "$RD/versao_iqtree.txt" 2>/dev/null || true)
[ -n "$v_concat" ] || v_concat=$v_cf
cat > "$S/versoes.yml" <<EOF
ASTER_WASTRAL:
  aster: 1.25.3.8
CLIPKIT:
  clipkit: 2.11.4
CONCATENA:
  python: 3.12.12
IQTREE_CONCAT:
  iqtree: ${v_concat}
IQTREE_CONCORDANCIA:
  iqtree: ${v_cf}
IQTREE_GENES:
  iqtree: 3.1.4
MAFFT_ALIGN:
  mafft: 7.520
  pigz: 2.8
PREPARA_LOCO:
  python: 3.12.12
EOF
cp -p "$S/versoes.yml" "$R/pipeline_info/versoes.yml"

docker_em "$S" "$IMG_PY" "python3 filogenia_util.py resumo \
  --params parametros.json --versoes versoes.yml --versao-nextflow 26.04.6 \
  --preparos preparos/* --alinhados alinhados/* --limpos limpos/* --relatorios-genes genes/*.iqtree \
  --dir-concat concat --dir-especies especies --cf cf/concordancia.gcf.cf.stat cf/concordancia.scfl.cf.stat \
  --saida-prefixo resumo"
cp -p "$S"/resumo.tsv "$S"/resumo.json "$S"/metodos.md "$R/"
echo "[pos] concluido ($(date '+%T'))"
