#!/usr/bin/env bash
# roda_fixo.sh — supermatrix 020210 rerun with the 48-group partition model that ModelFinder (MFP+MERGE) selected in the
# interrupted run. Resuming that run re-entered ModelFinder (IQ-TREE 3.1.4 re-tests single partitions whose cache entries
# hold only the merge-stage models, ~32 h), and a fixed-model resume from its checkpoint restarts the tree search anyway,
# so this is a clean run: same alignment, scheme, models, seed and options, minus the model search.
# Launched by lanca_fixo.ps1 as a scheduled task, so closing Claude Code does not SIGHUP it.
# $1 = IQ-TREE threads (default 8). A relaunch over an unfinished run resumes from concatenado.ckp.gz on its own: on
# 2026-10-19 the -T 8 run was stopped at iteration ~74 and relaunched with -T 18, after a test on a copy showed that the
# checkpoint, the 1000 UFBoot trees, the model parameters and the candidate set are all restored, not recomputed.
set -uo pipefail
NT=${1:-8}
AQUI=$(dirname "$(readlink -f "$0")")
TD=/home/manec/nf-work/filogenia/work/81/e3fdfb681f0da353c2d65ff89f587d
F=/home/manec/nf-work/filogenia/manual/020210/fixo
RD=$F/run
LOG=$F/supermatriz_fixa.log
IMG_IQ=quay.io/biocontainers/iqtree:3.1.4--h8471819_0

mkdir -p "$RD"
cd "$RD" || exit 1
[ -e supermatriz.fa ] || ln -s "$(readlink -f "$TD/supermatriz.fa")" supermatriz.fa
cp -p "$TD/concatenado.best_scheme.nex" "$TD/concatenado.best_scheme" .
# modelos_fixos.nex is generated before launch and checked against the original LnL/df; refuse to run without it.
[ -s modelos_fixos.nex ] || { echo "sem modelos_fixos.nex em $RD" >> "$LOG"; exit 1; }

echo "inicio: $(date '+%F %T') threads=$NT" >> "$LOG"
docker run --rm --name supermatriz_fixa_020210 -u "$(id -u):$(id -g)" -v /home/manec/nf-work:/home/manec/nf-work \
  -w "$RD" "$IMG_IQ" /bin/bash -c "
iqtree -version 2>&1 | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -n1 > versao_iqtree.txt
iqtree -p modelos_fixos.nex -B 1000 --bnni --alrt 1000 --abayes --seed 12345 -o Gallus_gallus -st AA \
  -s supermatriz.fa -pre concatenado -T $NT" >> "$LOG" 2>&1
st=$?
echo "iqtree fim: $(date '+%F %T') status=$st" >> "$LOG"
# With fixed models IQ-TREE writes no best_model.nex; the concordance step (-p) and the publish expect one.
[ -f concatenado.best_model.nex ] || cp -p modelos_fixos.nex concatenado.best_model.nex
if [ "$st" -eq 0 ] && [ -f concatenado.iqtree ]; then
  RD=$RD bash "$AQUI/pos_concat_020210.sh" >> "$LOG" 2>&1
  echo "pos_concat fim: $(date '+%F %T') status=$?" >> "$LOG"
fi
