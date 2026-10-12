#!/usr/bin/env bash
# bio-stack — clones de terceiros que nao sao pacote instalavel.
#
#   bash vendor/clonar.sh
#
# Idempotente. O conteudo de vendor/ NAO e versionado aqui (sao clones, e o
# .gitignore os exclui): este script e o artefato que reconstroi.
#
# BoneJ_Headless (harshamarupudi56, BSD-3-Clause) nao existe no PyPI nem no
# conda — e um conjunto de scripts que pilota o Fiji+BoneJ em modo headless
# pela API do ImageJ. O repositorio inteiro tem 415 MB, e quase tudo e imagem
# de exemplo em ROIs/ e Examples/; o codigo sao ~20 KB. Por isso o checkout e
# esparso: `--filter=blob:none` nao baixa blob que o checkout nao pedir, e
# sparse-checkout pede so os diretorios de codigo.
#
# Para trazer tambem os exemplos (os 400 MB), depois do clone:
#   git -C vendor/BoneJ_Headless sparse-checkout add Examples ROIs
set -euo pipefail

aqui="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# opencap-processing (stanfordnmbl, Apache-2.0) tambem nao e pacote: sao
# scripts que importam `opensim`. O ambiente que os roda e ~/envs/opencap no
# WSL, montado por bootstrap/06_wsl_envs_biomec.sh. Clone completo e raso: os
# 13 MB sao codigo e modelo de exemplo, nao ha o que podar.
if [ ! -d "$aqui/opencap-processing/.git" ]; then
  git clone --depth 1 https://github.com/stanfordnmbl/opencap-processing.git     "$aqui/opencap-processing"
else
  git -C "$aqui/opencap-processing" fetch --depth 1 origin
  git -C "$aqui/opencap-processing" reset --hard origin/HEAD
fi
echo "opencap-processing em $aqui/opencap-processing"
echo "  rodar com: pixi run --manifest-path ~/envs/opencap/pixi.toml python <script>"
echo "  o servico OpenCap exige conta: python createAuthenticationEnvFile.py"

if [ ! -d "$aqui/BoneJ_Headless/.git" ]; then
  git clone --filter=blob:none --no-checkout --depth 1 \
    https://github.com/harshamarupudi56/BoneJ_Headless.git \
    "$aqui/BoneJ_Headless"
  git -C "$aqui/BoneJ_Headless" sparse-checkout init --cone
  git -C "$aqui/BoneJ_Headless" sparse-checkout set \
    BoneJ_Scripts BoneJ_Analysis Secondary_Scripts Doc
  git -C "$aqui/BoneJ_Headless" checkout
else
  git -C "$aqui/BoneJ_Headless" fetch --depth 1 origin
  git -C "$aqui/BoneJ_Headless" reset --hard origin/HEAD
fi

echo "BoneJ_Headless em $aqui/BoneJ_Headless"
echo "Requer Fiji com a update site BoneJ habilitada. Nesta maquina:"
echo "  C:/Users/Manec/Apps/Fiji.app/ImageJ-win64.exe"
