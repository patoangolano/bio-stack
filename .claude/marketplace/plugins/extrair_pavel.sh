#!/usr/bin/env bash
# Materializa as 213 skills de Pavel-Kravchenko/Bioinformatics em bioinfo-213/.
#
# Por que archive e nao checkout: o indice do repositorio tem
# 'Legacy/Usefull_things_for_the_lab/<...>/Thumbs.db:encryptable'. O dois-pontos
# e sintaxe de alternate data stream no NTFS, e o git aborta a operacao INTEIRA
# com "invalid path" antes de olhar o filtro. Medido em 2026-10-11, nesta ordem:
#   git clone ............. arvore de trabalho vazia, exit 0, sem aviso
#   git checkout main ..... error: invalid path '<...>Thumbs.db:encryptable'
#   sparse-checkout set ... mesma falha: a validacao roda antes do filtro
#   git archive ........... mesma falha, "failed to unpack tree object"
#   archive com core.protectNTFS=false ... passa
# Como o archive leva so Skills/, o arquivo problematico nunca e materializado.
#
# Efeito colateral util: Course/ e Legacy/ ficam fora. O repositorio inteiro
# tem 496 MB; as 213 skills somam 2,5 MB — o peso era material de curso sem
# relacao com skill nenhuma.
set -euo pipefail
cd "$(dirname "$0")/bioinfo-213"

git fetch --quiet --depth 1 origin main
# Caminho relativo e fixo em vez de mktemp: sob Git Bash o mktemp devolve um
# caminho que o git e o tar traduzem de formas diferentes, e o tar reclama
# "This does not look like a tar archive" sobre um arquivo que existe.
tar_tmp=".pavel-archive.tar"
trap 'rm -f "$tar_tmp"' EXIT
git -c core.protectNTFS=false archive --format=tar -o "$tar_tmp" FETCH_HEAD Skills
tar -xf "$tar_tmp" -C .

# O repositorio chama de marketplace.json um arquivo que declara apenas
# {name, description, version, skills}. Sem 'owner' nem 'plugins' ele nao
# valida como marketplace, e presente no disco soma com a descoberta por
# convencao: 426 skills reportadas para 213 em disco.
rm -rf .claude-plugin

echo "skills extraidas: $(find Skills -name SKILL.md | wc -l)"
