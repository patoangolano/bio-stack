#!/usr/bin/env bash
# roda_topologia.sh — f07b_topologia.sh (constrained searches + KH/SH/AU) on the 020210 supermatrix, logging to
# topologia.log next to the hypotheses. Launched by lanca_topologia.ps1 as a scheduled task, so closing Claude Code
# does not SIGHUP it. f07b skips every hypothesis whose .treefile already exists: a relaunch picks up where it stopped.
# $1 = IQ-TREE threads (default 18).
set -uo pipefail
NT=${1:-18}
AQUI=$(dirname "$(readlink -f "$0")")
RUN=/mnt/c/Users/Manec/bio-stack/runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit
H=/mnt/c/Users/Manec/bio-stack/runs/rapinantes/filogenomica/topologia_AA
LOG=$H/topologia.log

echo "inicio: $(date '+%F %T') threads=$NT" >> "$LOG"
bash "$AQUI/f07b_topologia.sh" "$RUN" "$H" "$NT" >> "$LOG" 2>&1
st=$?
echo "fim: $(date '+%F %T') status=$st" >> "$LOG"
