#!/usr/bin/env bash
set -euo pipefail
SRC="/mnt/c/Users/user/bio-stack/runs/enamel_relics"
O="/mnt/c/Users/user/bio-stack/runs/enamel_relics/v101/stats"
TMP="/tmp/enamel_phaseC_$$"
mkdir -p "$TMP/enam" "$TMP/idp" "$TMP/out"
cp "$SRC/intervals_enamel"/interval_*.fna "$TMP/enam/"
cp "$SRC/intervals_idp"/idp_*.fna "$TMP/idp/"
cp "$O/queries/enamel_1000shuf.faa" "$TMP/enamel_q.faa"
cp "$O/queries/idp_1000shuf.faa" "$TMP/idp_q.faa"

run_one() {
  local dom="$1" seg="$2" sp="$3"
  local qfa subj outf
  if [[ "$dom" == "enamel" ]]; then
    qfa="$TMP/enamel_q.faa"
    subj="$TMP/enam/interval_${sp}.fna"
  else
    qfa="$TMP/idp_q.faa"
    subj="$TMP/idp/idp_${sp}.fna"
  fi
  outf="$TMP/out/${dom}_${seg}_${sp}.max.tsv"
  tblastn -query "$qfa" -subject "$subj" -evalue 10 -seg "$seg" -max_target_seqs 50 -num_threads 1           -outfmt '6 qseqid pident length evalue bitscore'     | awk -F'\t' '{
        q=$1; b=$5+0;
        if (!(q in best) || b > best[q]) { best[q]=b; pid[q]=$2; len[q]=$3; ev[q]=$4; }
      } END {
        for (q in best) printf "%s\t%s\t%s\t%s\t%.2f\n", q, pid[q], len[q], ev[q], best[q];
      }' > "$outf"
}
export -f run_one
export TMP

SPS=(Amis Ggan Cpor Acar Pbiv Pmur Cpic Cmyd Gevg Psin Ggal Apla Shab Tgut Dnov Scam Arow Tgua)
echo "[WSL] Starting 72 parallel tblastn jobs (18 sp x 2 domains x 2 SEG modes)..."
parallel -j 16 run_one {1} {2} {3} ::: enamel idp ::: yes no ::: "${SPS[@]}"
cp "$TMP/out"/*.max.tsv "$O/tblastn_raw/"
rm -rf "$TMP"
echo "[WSL] All 72 jobs completed."
