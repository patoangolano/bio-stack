#!/usr/bin/env bash
# Reconstroi os repositorios de skills de terceiros que o marketplace local
# declara em ../.claude-plugin/marketplace.json.
#
# Por que os clones moram AQUI e nao num vendor/ na raiz: o
# `claude plugin install` recusa uma entrada de marketplace cujo caminho saia do
# diretorio do marketplace. Medido em 2026-10-11, a mensagem e explicita e
# inclui juncao de diretorio na recusa — "an absolute, climbing,
# network-shaped, backslash-containing or link-traversing entry". Entao nem
# `source: ../../vendor/x` nem um `mklink /J` apontando para fora resolvem; a
# arvore precisa ser fisica e interna.
#
# Por que o marketplace local e nao `claude plugin marketplace add <repo>`:
#   - quatro dos cinco nao tem .claude-plugin/marketplace.json, e o add
#     responde "Marketplace file not found";
#   - o do Pavel tem um arquivo com esse nome, mas declara {name, description,
#     version, skills} — sem 'owner' nem 'plugins', nao valida como marketplace;
#   - o do choxos valida, mas declara componentes no marketplace.json E em cada
#     plugin.json com strict:false, e o Claude Code recusa carregar isso:
#     "conflicting manifests: both plugin.json and marketplace entry specify
#     components". Declarado aqui como ponteiro, sem componentes, ele carrega.
#
# As entradas deste marketplace NAO enumeram skills de proposito: enumerar e
# manter o diretorio skills/ faz as duas fontes somarem. O bioinfo-213 reportou
# 426 skills para 213 em disco antes de a enumeracao sair.
set -euo pipefail
cd "$(dirname "$0")"

# destino=repositorio
REPOS=(
  "bio-scientist-crispr=kexindon/claude-bio-scientist"
  "bioinfo-pipelines=YuliaNuzhnenko/bioinformatics-agent-skills"
  "biomedical-ai=zamushwani/biomedical-ai-skills"
  "arms-biostat=koaeraser/ARMS-Biostat"
  "biostat-agent=choxos/BiostatAgent"
)

for par in "${REPOS[@]}"; do
  destino="${par%%=*}"
  repo="${par#*=}"
  if [ -d "$destino/.git" ]; then
    git -C "$destino" pull --ff-only --quiet
    echo "atualizado  $destino"
  else
    git clone --quiet --depth 1 "https://github.com/$repo" "$destino"
    echo "clonado     $destino"
  fi
done

# O repositorio do choxos traz o proprio marketplace.json. Mantido no disco ele
# vira um marketplace aninhado dentro deste, com os componentes em conflito que
# o cabecalho descreve. Removido, os quatro plugin.json bastam.
rm -f biostat-agent/.claude-plugin/marketplace.json

# O do Pavel precisa de clone sem checkout: ver extrair_pavel.sh.
if [ ! -d bioinfo-213/.git ]; then
  git clone --quiet --depth 1 --no-checkout --filter=blob:limit=2m \
    https://github.com/Pavel-Kravchenko/Bioinformatics bioinfo-213
  echo "clonado     bioinfo-213 (sem checkout)"
fi
./extrair_pavel.sh

echo
for d in */; do
  d="${d%/}"
  printf '%-24s %6s  skills=%s\n' "$d" \
    "$(du -sh --exclude=.git "$d" | cut -f1)" \
    "$(find "$d" -name SKILL.md | wc -l)"
done
echo
echo "Agora: claude plugin marketplace update bio-stack"
echo "Depois, para cada plugin: claude plugin install <nome>@bio-stack --scope project"
