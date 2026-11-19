#!/usr/bin/env bash
# Phase 2.2 — Portable WSL2 bioinformatics alignment suite for the enamel-gene relics study.
# Runs:
#   1. Tool version logging + SEG low-complexity masking quantification (segmasker)
#   2. tBLASTn (-seg yes AND -seg no) of Alligator (full RefSeq + legacy) & Anolis queries
#      plus 20 composition-matched shuffles across all 18 NRG1-LPL enamel intervals
#   3. tBLASTn (-seg yes AND -seg no) of osteogenic IDP controls (IBSP, SPP1, SPARCL1)
#      plus 20 composition-matched shuffles across all 18 SIBLING/IDP intervals
#   4. miniprot protein-to-genome spliced alignment (ENAM full RefSeq, ENAM legacy, AMBN,
#      AMTN, and positive control LPL) across all 18 NRG1-LPL intervals
#   5. Whole-proteome BLAST database construction for Alligator mississippiensis (GCF_030867095.1)
#
# Usage (from Windows PowerShell):
#   wsl -d Ubuntu-26.04 -- bash -l /mnt/c/Users/user/bio-stack/scripts/enamel_relics/02_run_alignments_wsl.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO_ROOT/runs/enamel_relics"
QDIR="$OUT/queries"
ENAM_DIR="$OUT/intervals_enamel"
IDP_DIR="$OUT/intervals_idp"
REF_DIR="$OUT/refdata"

mkdir -p "$OUT/tblastn_enamel_seg_yes" \
         "$OUT/tblastn_enamel_seg_no" \
         "$OUT/tblastn_idp_seg_yes" \
         "$OUT/tblastn_idp_seg_no" \
         "$OUT/miniprot" \
         "$OUT/blastdb" \
         "$OUT/segmask"

THREADS=8
FMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

echo "=== [1/5] Recording exact tool versions ==="
{
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "uname=$(uname -a)"
  echo "tblastn=$(tblastn -version | head -n 1)"
  echo "blastx=$(blastx -version | head -n 1)"
  echo "blastn=$(blastn -version | head -n 1)"
  echo "segmasker=$(segmasker -version 2>&1 | head -n 1 || echo unavailable)"
  echo "diamond=$(diamond version 2>&1 | head -n 1 || echo unavailable)"
  echo "miniprot=$(miniprot --version 2>&1 | head -n 1)"
  echo "mafft=$(mafft --version 2>&1 | head -n 1)"
  echo "muscle=$(muscle -version 2>&1 | head -n 1)"
  echo "samtools=$(samtools --version | head -n 1)"
  echo "seqkit=$(seqkit version 2>&1 | head -n 1)"
  echo "datasets=$(datasets --version 2>&1 | head -n 1)"
} > "$OUT/tool_versions.txt"
cat "$OUT/tool_versions.txt"

echo ""
echo "=== [1b] Quantifying SEG low-complexity masking on real query proteins ==="
# Extract non-shuffled IDP queries into idp_real.faa
awk '/^>/{p=($0 !~ /~shuf/)} p' "$QDIR/idp_with_shuffles.faa" > "$QDIR/idp_real.faa"
cat "$QDIR/enamel_real.faa" "$QDIR/idp_real.faa" > "$QDIR/all_real_queries.faa"
if command -v segmasker >/dev/null 2>&1; then
  segmasker -in "$QDIR/all_real_queries.faa" -infmt fasta -outfmt interval \
            -window 12 -locut 2.2 -hicut 2.5 \
            -out "$OUT/segmask/seg_intervals.txt"
  echo "  segmasker intervals written to $OUT/segmask/seg_intervals.txt"
else
  echo "  segmasker not found in PATH"
fi

echo ""
echo "=== [2/5] Running tBLASTn (-seg yes and -seg no) on 18 NRG1-LPL enamel intervals ==="
for f in "$ENAM_DIR"/interval_*.fna; do
  sp="$(basename "$f" .fna)"
  sp="${sp#interval_}"
  tblastn -query "$QDIR/enamel_with_shuffles.faa" -subject "$f" \
          -evalue 10 -seg yes -max_target_seqs 50 -num_threads "$THREADS" \
          -outfmt "$FMT" > "$OUT/tblastn_enamel_seg_yes/${sp}.tsv"
  tblastn -query "$QDIR/enamel_with_shuffles.faa" -subject "$f" \
          -evalue 10 -seg no -max_target_seqs 50 -num_threads "$THREADS" \
          -outfmt "$FMT" > "$OUT/tblastn_enamel_seg_no/${sp}.tsv"
  printf '  [ENAMEL] %-4s  seg_yes=%5d HSPs  seg_no=%5d HSPs\n' \
         "$sp" "$(wc -l < "$OUT/tblastn_enamel_seg_yes/${sp}.tsv")" \
         "$(wc -l < "$OUT/tblastn_enamel_seg_no/${sp}.tsv")"
done

echo ""
echo "=== [3/5] Running tBLASTn (-seg yes and -seg no) on 18 osteogenic IDP intervals ==="
for f in "$IDP_DIR"/idp_*.fna; do
  sp="$(basename "$f" .fna)"
  sp="${sp#idp_}"
  tblastn -query "$QDIR/idp_with_shuffles.faa" -subject "$f" \
          -evalue 10 -seg yes -max_target_seqs 50 -num_threads "$THREADS" \
          -outfmt "$FMT" > "$OUT/tblastn_idp_seg_yes/${sp}.tsv"
  tblastn -query "$QDIR/idp_with_shuffles.faa" -subject "$f" \
          -evalue 10 -seg no -max_target_seqs 50 -num_threads "$THREADS" \
          -outfmt "$FMT" > "$OUT/tblastn_idp_seg_no/${sp}.tsv"
  printf '  [IDP   ] %-4s  seg_yes=%5d HSPs  seg_no=%5d HSPs\n' \
         "$sp" "$(wc -l < "$OUT/tblastn_idp_seg_yes/${sp}.tsv")" \
         "$(wc -l < "$OUT/tblastn_idp_seg_no/${sp}.tsv")"
done

echo ""
echo "=== [4/5] Running miniprot across all 18 NRG1-LPL intervals ==="
# Parameters match previous calibrated run: -j 2 -k 4 -M 0 -n 1 -l 4 --outc=0.01 --outs=0.5
MP_OPTS=(-j 2 -k 4 -M 0 -n 1 -l 4 --outc=0.01 --outs=0.5 -t "$THREADS")
for g in ENAM_Amis ENAM_Amis_legacy AMBN_Amis AMTN_Amis AMTN_Amis_gb LPL_Amis; do
  q="$QDIR/mp_${g}.faa"
  for f in "$ENAM_DIR"/interval_*.fna; do
    sp="$(basename "$f" .fna)"
    sp="${sp#interval_}"
    miniprot "${MP_OPTS[@]}" "$f" "$q" 2>/dev/null > "$OUT/miniprot/mp_${g}_${sp}.paf" || true
  done
  n_hit=$(find "$OUT/miniprot" -name "mp_${g}_*.paf" -size +0c | wc -l)
  printf '  [MINIPROT] %-16s  aligned in %2d / 18 species\n' "$g" "$n_hit"
done

echo ""
echo "=== [5/5] Building whole-proteome BLAST database for Alligator (GCF_030867095.1) ==="
makeblastdb -in "$REF_DIR/Amis_protein.faa" -dbtype prot \
            -out "$OUT/blastdb/amis_proteome" >/dev/null
n_prot=$(grep -c '^>' "$REF_DIR/Amis_protein.faa")
echo "  Alligator proteome BLAST DB ready: $n_prot proteins"
echo "=== Phase 2.2 WSL alignments complete ==="
