#!/usr/bin/env bash
# exportar_repositorio.sh — mirrors the rapinantes part of bio-stack into the project's own repository
# (github.com/patoangolano/filogenia-rapinantes). bio-stack stays the source: this copies the files tracked at one
# bio-stack commit, replaces everything in the destination clone and commits there, citing the source commit.
# Uso: bash scripts/rapinantes/exportar_repositorio.sh [destino=../filogenia-rapinantes] [commit=HEAD]
set -euo pipefail
BIO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
DEST=${1:-$BIO/../filogenia-rapinantes}
REV=$(git -C "$BIO" rev-parse "${2:-HEAD}")
CAMINHOS=(docs/metodos_rapinantes.md reports research_notes scripts/rapinantes runs/rapinantes
          runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit
          runs/filogenia/20261018-084205-rapinantes_DNA12
          runs/filogenia/20261018-112050-rapinantes_AA_famsa_clipkit)
[ -d "$DEST/.git" ] || { echo "sem repositorio git em $DEST" >&2; exit 1; }

# Mirror, not merge: a file deleted in bio-stack must go here too. Only tracked files travel, so whatever bio-stack's
# .gitignore keeps out (per-gene figures, third-party downloads, checkpoints) stays out of this repository as well.
git -C "$DEST" rm -rq --ignore-unmatch .
git -C "$BIO" archive "$REV" -- "${CAMINHOS[@]}" | tar -x -C "$DEST"
git -C "$BIO" show "$REV:scripts/rapinantes/repositorio/README.md" > "$DEST/README.md"
printf '* text=auto eol=lf\n*.ps1 text eol=crlf\n' > "$DEST/.gitattributes"

git -C "$DEST" add -A
if git -C "$DEST" diff --cached --quiet; then
  echo "nada mudou desde a ultima exportacao"
  exit 0
fi
git -C "$DEST" commit -q -m "sincroniza com bio-stack ${REV:0:7}" \
  -m "Fonte: https://github.com/patoangolano/bio-stack/tree/$REV"
git -C "$DEST" log --oneline -n 1
